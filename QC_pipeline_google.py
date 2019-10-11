
'''
QC for 1 chromosome
author: Pavlos Antoniou
date: 22/07/19
'''

import os
import hail as hl

import json

CHROMOSOME="chrY"

project_root = os.path.dirname(os.path.dirname(__file__))
print(project_root)


storage = os.path.join(project_root , "config_files/storage.json")

thresholds = os.path.join(project_root, "config_files/thresholds.json")

with open(f"{storage}", 'r') as f:
    storage = json.load(f)

with open(f"{thresholds}", 'r') as f:
    thresholds = json.load(f)

BUCKET = storage['intervalwgs']['google']['bucket']
#Define chromosome here
CHROMOSOME="chr1"
tmp_dir=storage['intervalwgs']['google']['tmp_dir']

if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    #Define the hail persistent storage directory
    hl.init(default_reference="GRCh38", tmp_dir=tmp_dir)

    VQSLOD_snps = hl.import_table(f"{BUCKET}/qc-files/VQSLOD_snps.bgz",
                                  types={"Locus": "locus<GRCh38>", "VQSLOD": hl.tfloat64})
    VQSLOD_indels = hl.import_table(f"{BUCKET}/qc-files/VQSLOD_indels.bgz",
                                    types={"Locus": "locus<GRCh38>", "VQSLOD": hl.tfloat64})
    sample_QC_nonHail = hl.import_table(f"{BUCKET}/qc-files/INTERVAL_WGS_Sample_QC_28-05-2019_fixed_header.txt",
                                        impute=True)
    gws_gwa_map = hl.import_table(f"{BUCKET}/qc-files/WGS-2-GWA_omicsMap.txt", impute=True)

    full_blood_count = hl.import_table(f"{BUCKET}/qc-files/Int50kSampleId.tsv", impute=True, key='Wgs_RAW_bl')


    print('Joining annotations')
    ja = sample_QC_nonHail.key_by('ID')
    gws_gwa_map = gws_gwa_map.key_by('unique_WGS_ID')

    ja = ja.annotate(gws_gwa_map=gws_gwa_map[ja['ID']])
    fbc = full_blood_count.key_by('Wgs_RAW_bl')
    ja = ja.annotate(fbc=full_blood_count[ja.gws_gwa_map.Wgs_RAW_bl])


    print(f'Read {CHROMOSOME} mt ')
    mt = hl.read_matrix_table(f"{BUCKET}/{CHROMOSOME}.mt")
    ####################################################################
    print('Split multialleles')
    mt_split1=hl.split_multi_hts(mt, keep_star=False)
    print('checkpoint split matrixtable')
    mt_split = mt_split1.checkpoint( f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-split-multi.mt", overwrite= True)

    #print('coalesce partition')
    # #mt_split = mt_split1.naive_coalesce(10000)
    #print('coalesce checkpoint')
    #mt_split= mt_split.checkpoint(f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}_coalesce_checkpoint.mt", overwrite = True)

    print('Annotate matrixtable with snp indel and other information')
    mt = mt_split.annotate_rows(Variant_Type=hl.cond((hl.is_snp(mt_split.alleles[0], mt_split.alleles[1])), "SNP",
                                                         hl.cond(
                                                             hl.is_insertion(mt_split.alleles[0], mt_split.alleles[1]),
                                                             "INDEL",
                                                             hl.cond(hl.is_deletion(mt_split.alleles[0],
                                                                                    mt_split.alleles[1]), "INDEL",
                                                                     "Other"))))

    #Unfiltered data summary stats:
    print('Unfiltered sample qc and write to google bucket')
    mt_sqc1_unfiltered = mt.annotate_cols(sample_QC_nonHail=sample_QC_nonHail.key_by("ID")[mt.s])
    mt_sqc2_unfiltered = hl.sample_qc(mt_sqc1_unfiltered, name='sample_QC_Hail')
    panda_df_unfiltered_table = mt_sqc2_unfiltered.cols().flatten()
    print('Export to table')
    panda_df_unfiltered_table.export(f"{BUCKET}/output-tables/{CHROMOSOME}/{CHROMOSOME}-sampleQC_unfiltered.tsv.bgz", header=True)

    print('VQSLOD scores filtering')
    mt = mt.annotate_rows(VQSLOD_SNP=VQSLOD_snps.key_by("Locus")[mt.locus])
    mt = mt.annotate_rows(VQSLOD_INDEL=VQSLOD_indels.key_by("Locus")[mt.locus])
    vqslod_filtered = (
            hl.case()
                .when((mt.Variant_Type == "SNP"), (mt.VQSLOD_SNP.VQSLOD >= thresholds['intervalwgs']['snp_vqsr_threshold']))
                .when((mt.Variant_Type == "INDEL"), (mt.VQSLOD_INDEL.VQSLOD >= thresholds['intervalwgs']['indel_vqsr_threshold']))
                .default(False)  # remove everything else
        )

    mt_vqslod_filtered = mt.filter_rows(vqslod_filtered)

    print('Vqslod checkpoint')
    mt1 = mt_vqslod_filtered.checkpoint(
            f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}_vqslod_filtered_checkpoint.mt",  overwrite= True)




    print('Sample QC filtering ')
    mt_sqc1 = mt1.annotate_cols(
        sample_QC_nonHail=sample_QC_nonHail.key_by("ID")[mt1.s])

    mt_sqc1_filtered = mt_sqc1.filter_cols(
        (mt_sqc1.sample_QC_nonHail.PASS_Depth == 1) &
        (mt_sqc1.sample_QC_nonHail.PASS_ID == 1) &
        (mt_sqc1.sample_QC_nonHail.PASS_Median_FreeMix == 1) &
        (mt_sqc1.sample_QC_nonHail.PASS_NRD == 1) &
        (mt_sqc1.sample_QC_nonHail.PASS_SampleSwap == 1) &
        (mt_sqc1.sample_QC_nonHail.PASS_Sex == 1) &
        (mt_sqc1.sample_QC_nonHail.PASS_DUP == 1)
    )

    mt_sqc2 = hl.sample_qc(mt_sqc1_filtered, name='sample_QC_Hail')
    print('write to to google bucket')
    mt_sqc2 = mt_sqc2.checkpoint(f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-sampleQC_filtered.mt", overwrite= True)
    filter_sampleqc_table=mt_sqc2.cols().flatten()
    filter_sampleqc_table.export(f"{BUCKET}/output-tables/{CHROMOSOME}/{CHROMOSOME}-sampleQC_filtered.tsv.bgz", header=True)

##########################################
    #Variant QC
    print('Variant QC and write matrixtable to to google bucket:')
    mt_sqc_vqc = hl.variant_qc(mt_sqc2, name='variant_QC_Hail')

    mt_sqc_vqc_filtered = mt_sqc_vqc.filter_rows(
        (mt_sqc_vqc.variant_QC_Hail.call_rate >= 0.98) &
        (mt_sqc_vqc.variant_QC_Hail.p_value_hwe >= 10 ** -6))


    # drop those structure run additional sample qc and variant qc and have the final values
    # sample and filter qc tab tables
    fields_to_drop = ['variant_QC_Hail', 'sample_QC_Hail']

    mt1 = mt_sqc_vqc_filtered.drop(*fields_to_drop)

    mt2 = hl.sample_qc(mt1,name='sample_QC_Hail')
    mt3 = hl.variant_qc(mt2, name='variant_QC_Hail')

    mt3 = mt3.checkpoint(
            f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-sampleqc-variantqc-filtered-FINAL.mt", overwrite= True)

    mt3_rows = mt3.rows()
    mt3_rows.select(mt3_rows.variant_QC_Hail).flatten().export(f"{BUCKET}/output-tables/{CHROMOSOME}/{CHROMOSOME}-variantQC-sampleQC_filtered_FINAL.tsv.bgz",
                                                               header=True)
