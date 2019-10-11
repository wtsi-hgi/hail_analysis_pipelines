
'''
QC for 1 chromosome
author: Pavlos Antoniou
date: 22/07/19
'''

import os
import hail as hl
import pyspark
import json

CHROMOSOME="chrY"

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
    VQSLOD_snps = hl.import_table(storage["intervalwgs"]["s3"]["vqsr_snp"],
                                  types={"Locus": "locus<GRCh38>", "VQSLOD": hl.tfloat64})
    VQSLOD_indels = hl.import_table(storage["intervalwgs"]["s3"]["vqsr_indel"],
                                    types={"Locus": "locus<GRCh38>", "VQSLOD": hl.tfloat64})
    sample_QC_nonHail = hl.import_table(storage.intervalwgs.s3.sampleQC_non_hail, impute=True)


    print(f"Reading {CHROMOSOME} mt")
    mt = hl.read_matrix_table(f"{storage['s3']['matrixtables']}/{CHROMOSOME}.mt")

    print("Splitting mt and writing out split mt")
    mt_split = hl.split_multi_hts(mt, keep_star=False)


    mt_split = mt_split.checkpoint(f"{tmp_dir}/matrixtable/{CHROMOSOME}-split-multi.mt",  overwrite=True)
    print("Finished splitting and writing mt. ")

    #print("Repartitioning to 1000 partitions ")
    #mt_split = mt_split.naive_coalesce(1000)





    print('Annotating rows with snp and indel info')
    mt = mt_split.annotate_rows(
                Variant_Type=hl.cond((hl.is_snp(mt_split.alleles[0], mt_split.alleles[1])), "SNP",
                                     hl.cond(
                                         hl.is_insertion(mt_split.alleles[0], mt_split.alleles[1]),
                                         "INDEL",
                                         hl.cond(hl.is_deletion(mt_split.alleles[0],
                                                                mt_split.alleles[1]), "INDEL",
                                                 "Other"))))

            # Unfiltered data summary stats:
    print("Finished annotating rows, annotating columns now")
    mt_sqc1_unfiltered = mt.annotate_cols(sample_QC_nonHail=sample_QC_nonHail.key_by("ID")[mt.s])
    mt_sqc2_unfiltered = hl.sample_qc(mt_sqc1_unfiltered, name='sample_QC_Hail')
    panda_df_unfiltered_table = mt_sqc2_unfiltered.cols().flatten()

    print("Outputting table of sample qc")
    panda_df_unfiltered_table.export(f"{tmp_dir}/output-tables/{CHROMOSOME}_sampleQC_unfiltered.tsv.bgz", header=True)

    mt2 = hl.variant_qc(mt, name='variant_QC_Hail')

    print('Exporting variant qc pandas table to disk')
    mt_rows = mt2.rows()
    mt_rows.select(mt_rows.variant_QC_Hail).flatten().export(f"{tmp_dir}/output-tables/{CHROMOSOME}-variantQC_unfiltered_pandaDF.tsv.bgz",
                                                             header=True)

    print("Annotation VQSLOD snp and indel scores")
    mt = mt.annotate_rows(VQSLOD_SNP=VQSLOD_snps.key_by("Locus")[mt.locus])
    mt = mt.annotate_rows(VQSLOD_INDEL=VQSLOD_indels.key_by("Locus")[mt.locus])

    print("Filtering on VQSLOD scores")
    vqslod_filtered = (
                hl.case()
                    .when((mt.Variant_Type == "SNP"), (mt.VQSLOD_SNP.VQSLOD >= thresholds['intervalwgs']['snp_vqsr_threshold']))
                    .when((mt.Variant_Type == "INDEL"), (mt.VQSLOD_INDEL.VQSLOD >= thresholds['intervalwgs']['indel_vqsr_threshold']))
                    .default(False)  # remove everything else
            )

    mt_vqslod_filtered = mt.filter_rows(vqslod_filtered)
    print("Finished filtering. Writing out matrixtable...")

    mt1 = mt_vqslod_filtered.checkpoint(f"{tmp_dir}/checkpoints/{CHROMOSOME}_vqslod_filtered_checkpoint.mt", overwrite=True)
    print("Finished writing vqslod filtered matrixtable")

    print("Annotating columns with sample qc")
    mt_sqc1 = mt1.annotate_cols(
                sample_QC_nonHail=sample_QC_nonHail.key_by("ID")[mt1.s])
    print("Filtering on sample qc")
    mt_sqc1_filtered = mt_sqc1.filter_cols(
                (mt_sqc1.sample_QC_nonHail.PASS_Depth == 1 ) &
                (mt_sqc1.sample_QC_nonHail.PASS_ID == 1) &
                (mt_sqc1.sample_QC_nonHail.PASS_Median_FreeMix == 1) &
                (mt_sqc1.sample_QC_nonHail.PASS_NRD == 1) &
                (mt_sqc1.sample_QC_nonHail.PASS_SampleSwap == 1) &
                (mt_sqc1.sample_QC_nonHail.PASS_Sex == 1) &
                (mt_sqc1.sample_QC_nonHail.PASS_DUP == 1)
            )
    print("Writing out filtered sample qc checkpoint")
    mt_sqc2 = hl.sample_qc(mt_sqc1_filtered, name='sample_QC_Hail')

    mt_sqc2 = mt_sqc2.checkpoint(f"{tmp_dir}/checkpoints/{CHROMOSOME}_sampleQC_step1_filtered_checkpoint.mt",
                                          overwrite=True)
    filter_sampleqc_table = mt_sqc2.cols().flatten()
    filter_sampleqc_table.export(f"{tmp_dir}/output-tables/{CHROMOSOME}-sampleQC_filtered.tsv.bgz", header=True)


    print("Variant qc:")
    mt_sqc_vqc = hl.variant_qc(mt_sqc2, name='variant_QC_Hail')
    mt_sqc_vqc_filtered = mt_sqc_vqc.filter_rows(
                (mt_sqc_vqc.variant_QC_Hail.call_rate >= 0.98) &
                (mt_sqc_vqc.variant_QC_Hail.p_value_hwe >= 10 ** -6))


    fields_to_drop = ['variant_QC_Hail', 'sample_QC_Hail']

    mt1 = mt_sqc_vqc_filtered.drop(*fields_to_drop)

    mt2 = hl.sample_qc(mt1, name='sample_QC_Hail')
    mt3 = hl.variant_qc(mt2, name='variant_QC_Hail')

    mt3 = mt_sqc_vqc_filtered.checkpoint(
        f"{tmp_dir}/checkpoints/{CHROMOSOME}-full-sampleqc-variantqc-filtered-FINAL.mt", overwrite=True)
    mt3_rows = mt3.rows()
    mt3_rows.select(mt3_rows.variant_QC_Hail).flatten().export(
        f"{tmp_dir}/output-tables/{CHROMOSOME}-variantQC_sampleQC_filtered_FINAL.tsv.bgz", header=True)



