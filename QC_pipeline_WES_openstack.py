
'''
QC for WES project
author: Pavlos Antoniou
date: 11/10/2019
'''

import os
import hail as hl
import pyspark
import json



project_root = os.path.dirname(os.path.dirname(__file__))
print(project_root)

s3credentials = os.path.join(project_root, "config_files/s3_credentials.json")
print(s3credentials)

storage = os.path.join(project_root , "config_files/storage.json")

thresholds = os.path.join(project_root, "config_files/thresholds.json")

with open(f"{s3credentials}", 'r') as f:
    credentials = json.load(f)

with open(f"{storage}", 'r') as f:
    storage = json.load(f)

with open(f"{thresholds}", 'r') as f:
    thresholds = json.load(f)


if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    #Define the hail persistent storage directory
    tmp_dir = os.path.join(os.environ["HAIL_HOME"], "tmp")
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    #s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])

    #1. Import required files
    print("1. import required tables in hail for the project.")

    exclude_samples_table = hl.import_table(storage["intervalwes"]["s3"]["exclude_samples"], no_header=True, key='f0')
    input_vcf=storage["intervalwes"]["s3"]["exome"]
    VQSLOD_snps = hl.import_table(storage["intervalwes"]["s3"]["vqsr_snp"],
                                  types={"Locus": "locus<GRCh38>", "VQSLOD": hl.tfloat64})
    VQSLOD_indels = hl.import_table(storage["intervalwes"]["s3"]["vqsr_indel"],
                                    types={"Locus": "locus<GRCh38>", "VQSLOD": hl.tfloat64})
    #1. Import VCF/ MT
    #print("1. Import VCF")
    #mt = hl.import_vcf(input_vcf,force_bgz=True, reference_genome='GRCh38', skip_invalid_loci=True)
    mt = hl.read_matrix_table(f"{tmp_dir}/matrixtables/WES.mt")
    #2. Remove samples that have not passed initial QC:
    print("2. Remove samples that have not passed initial QC:")
    mt_result = mt.filter_cols(hl.is_defined(exclude_samples_table[mt.s]), keep=False)

    #3. Split multi
    print("3. Split multi")
    mt_split = hl.split_multi_hts(mt_result, keep_star=False)

    # 4. annotate SNPs,indels
    print('Annotating rows with snp and indel info')
    mt = mt_split.annotate_rows(
        Variant_Type=hl.cond((hl.is_snp(mt_split.alleles[0], mt_split.alleles[1])), "SNP",
                             hl.cond(
                                 hl.is_insertion(mt_split.alleles[0], mt_split.alleles[1]),
                                 "INDEL",
                                 hl.cond(hl.is_deletion(mt_split.alleles[0],
                                                        mt_split.alleles[1]), "INDEL",
                                         "Other"))))


    #4. Sample qc and variant qc
    print("4. Sample qc and variant qc ")
    mt_sampleqc = hl.sample_qc(mt, name='sample_QC_Hail')
    mt2 = hl.variant_qc(mt_sampleqc, name='variant_QC_Hail')


    #5.Annotate COMMON AND RARE VARIANTS to apply separate filters
    print("Annotate COMMON AND RARE VARIANTS to apply separate filters")
    #mt_common = mt_filtered.filter_rows(mt_filtered.variant_qc.AF[1] > 0.05)
    mt2 = mt2.annotate_rows(maf=hl.cond(mt2.variant_QC_Hail.AF[1] < 0.01, "< 1%",
                                        hl.cond(mt2.variant_QC_Hail.AF[1] < 0.05, "1%-5%", ">5%")))


    #6. Common variants  filtering:
    print("6. Common variants  filtering:")
    mt=mt2
    mt_filtered_variants_common = mt.filter_rows(
        (mt.maf == "< 1%") | #let all rare variants pass
        (
        (mt.maf != "< 1%") &
        ((mt.variant_QC_Hail.p_value_hwe > 10 ** -5) & (mt.variant_QC_Hail.call_rate >= 0.97))
        )
    )


    #VQSR

    mt=mt_filtered_variants_common
    ######## VQSR filtering

    print("Annotation VQSLOD snp and indel scores")
    mt = mt.annotate_rows(VQSLOD_SNP=VQSLOD_snps.key_by("Locus")[mt.locus])
    mt = mt.annotate_rows(VQSLOD_INDEL=VQSLOD_indels.key_by("Locus")[mt.locus])

    print("Filtering on VQSLOD scores")
    vqslod_filtered = (
        hl.case()
            .when((mt.Variant_Type == "SNP"), (mt.VQSLOD_SNP.VQSLOD >= thresholds['intervalwes']["snps"]['snp_vqsr_threshold']))
            .when((mt.Variant_Type == "INDEL"),
                  (mt.VQSLOD_INDEL.VQSLOD >= thresholds['intervalwes']["indels"]['indel_vqsr_threshold']))
            .default(False)  # remove everything else
    )

    mt_vqslod_filtered = mt.filter_rows(vqslod_filtered)
    print("Finished vqslod filtering. Writing out matrixtable...")

    mt_vqslod_filtered = mt_vqslod_filtered.checkpoint(
        f"{tmp_dir}/intervalwes/WES_vqslod_filtered_checkpoint.mt", overwrite=True)
    print("Finished writing vqslod filtered matrixtable")



    #7.Rare variants filtering
    print("Rare variants filtering")

    mt = mt_vqslod_filtered


    mt_entries_filtered = mt.filter_entries(
        (mt.maf != "< 1%") |
        (
                (mt.Variant_Type == "Other") |
                ((mt.Variant_Type == "SNP") & (mt.GQ >= thresholds["intervalwes"]["snps"]["GQ"]) & (mt.DP >= thresholds["intervalwes"]["snps"]["DP"])) |
                ((mt.Variant_Type == "INDEL") & (mt.GQ >= thresholds["intervalwes"]["indels"]["GQ"]) & (mt.DP >= thresholds["intervalwes"]["snps"]["DP"] ))
        )
    )

    ###############Allelic Balance

    mt_sqc_vqc = mt_entries_filtered.filter_entries(hl.is_defined(mt_entries_filtered.GT))

    ab = mt_sqc_vqc.AD[1] / hl.sum(mt_sqc_vqc.AD)

    filter_condition_ab = (
        hl.case(missing_false=True)
            #     .when(((mt_sqc_vqc.GT.is_hom_ref()) & (mt_sqc_vqc.working_freq == "Rare")), ab > 0.1)
            .when(((mt_sqc_vqc.GT.is_het()) & (mt_sqc_vqc.maf != "< 1%")), (ab < 0.15) | (ab > 0.80))
            #     .when(((mt_sqc_vqc.GT.is_hom_var()) & (mt_sqc_vqc.working_freq == "Rare")), ab < 0.9)
            .default(False)  # remove everything else
    )

    mt_sqc_vqc_GT = mt_sqc_vqc.filter_entries(filter_condition_ab, keep=False)

    #Now need to filter rows again by doing a new variant_qc.
    #Drop the previous variant_qc and sample_qc
    print("8. Drop old sample and variant QC and calculate new ones after filtering")

    fields_to_drop = ['variant_QC_Hail', 'sample_QC_Hail']

    mt1 = mt_sqc_vqc_GT.drop(*fields_to_drop)

    mt2 = hl.sample_qc(mt1, name='sample_QC_Hail')
    mt3 = hl.variant_qc(mt2, name='variant_QC_Hail')

    #Remove variants with missingness > 3%
    mt= mt3.filter_rows(mt3.variant_QC_Hail.call_rate >= 0.97)

    print("Finished filtering. Now writing out.")

    #8.Export VCF only variants -no genotypes
    print("Export VCF only variants -no genotypes")
    mt1 = mt.select_entries()
    mt2 = mt1.filter_cols(mt1['s'] == 'samplenone')
    hl.export_vcf(mt2, f"{tmp_dir}/intervalwes/VCFs/exome.vcf.bgz")

    #9. Write matrixtable
    print("Write matrixtable")
    mt = mt.checkpoint(f"{tmp_dir}/intervalwes/exome_filtered.mt",  overwrite=True)
    print("Finished writing mt. ")




