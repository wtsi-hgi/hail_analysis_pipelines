
'''
QC for 1 chromosome
author: Pavlos Antoniou
date: 22/07/19
'''

import os
import hail as hl
import pyspark
import json
import sys




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


    CHROMOSOME="chr20"

    mt_sqc2 = hl.read_matrix_table(f"{tmp_dir}/checkpoints/{CHROMOSOME}_sampleQC_filtered_checkpoint.mt")

    ##### ADD allele bias script here

    initial_geno = mt_sqc2.aggregate_entries(hl.agg.fraction(hl.is_defined(mt_sqc2.GT)))

    print(f'Defined genotypes: {initial_geno * 100:.2f}%.')
    initial_missing = 100 - (initial_geno * 100)
    print(f'Initial missing genotype: {initial_missing:.2f}%.')

    ## Remove initial missing genotypes
    mt_sqc2 = mt_sqc2.filter_entries(hl.is_defined(mt_sqc2.GT))

    ab = mt_sqc2.AD[1] / hl.sum(mt_sqc2.AD)

    filter_condition_ab = (
        hl.case(missing_false=True)
            .when(mt_sqc2.GT.is_hom_ref(), ab > 0.1)
            .when(mt_sqc2.GT.is_het(), (ab < 0.20) | (ab > 0.80))
            .when(mt_sqc2.GT.is_hom_var(), ab < 0.9)
            .default(False)  # remove everything else
    )

    fraction_filtered = mt_sqc2.aggregate_entries(hl.agg.fraction(filter_condition_ab))
    print(f'Filtering {fraction_filtered * 100:.2f}% entries out of downstream analysis.')

    print(
        f'Total Filtering {(fraction_filtered * 100) + initial_missing:.2f}% entries including initial missing data out of downstream analysis.')

    Total_geno = mt_sqc2.aggregate_entries(hl.agg.count_where(hl.is_defined(mt_sqc2.GT)))
    print('Total genotypes: ' + str(Total_geno))

    mt_sqc2_GT = mt_sqc2.filter_entries(filter_condition_ab, keep=False)


    ############################

    pro_AD_DP = hl.sum(mt_sqc2_GT.AD) / mt_sqc2_GT.DP

    Other_filters = (
        hl.case(missing_false=True)
            .when(hl.is_defined(mt_sqc2_GT.GT), (mt_sqc2_GT.DP > 100) | (pro_AD_DP < 0.9))
            .default(False)  # remove everything else
    )

    mt_sqc3 = mt_sqc2_GT.filter_entries(Other_filters, keep=False)

    ## Saving on s3
    #mt_sqc2 = mt_sqc2_GT.checkpoint('s3a://interval-wgs/checkpoints_new/chr20_sampleQC_step1_filtered_allele_balance_checkpoint.mt', _read_if_exists = True)

    ## Saving on local volume
    mt_sqc3= mt_sqc3.checkpoint(f"{tmp_dir}/checkpoints/{CHROMOSOME}_sampleQC_step1_filtered_GenotypeQC_checkpoint.mt",
                                    overwrite=True)

    ######################
    ##### VARIANT qc
    ######################
    print("Variant qc:")
    mt_sqc_vqc = hl.variant_qc(mt_sqc3, name='variant_QC_Hail')
    mt_sqc_vqc_filtered = mt_sqc_vqc.filter_rows(
        (mt_sqc_vqc.variant_QC_Hail.call_rate >= 0.95) &
        (mt_sqc_vqc.variant_QC_Hail.p_value_hwe >= 10 ** -6) &
        (mt_sqc_vqc.variant_QC_Hail.gq_stats.mean >= 20) &
        (mt_sqc_vqc.variant_QC_Hail.AC[1] >= 1)
    )

    #####################################################################
    ###################### FINAL QC AFTER FILTERING  ####################
    #####################################################################

    fields_to_drop = ['variant_QC_Hail', 'sample_QC_Hail']

    mt1 = mt_sqc_vqc_filtered.drop(*fields_to_drop)



    mt2 = hl.sample_qc(mt1, name='sample_QC_Hail')
    mt3 = hl.variant_qc(mt2, name='variant_QC_Hail')

    mt3 = mt3.checkpoint(
        f"{tmp_dir}/checkpoints/{CHROMOSOME}-full-sampleqc-variantqc-filtered-FINAL.mt", overwrite=True)

    mt3_cols = mt3.cols()
    mt3_cols.flatten().export(
        f"{tmp_dir}/output-tables/{CHROMOSOME}-sampleQC_filtered_FINAL.tsv.bgz", header=True)

    mt3_rows = mt3.rows()
    mt3_rows.select(mt3_rows.variant_QC_Hail).flatten().export(
        f"{tmp_dir}/output-tables/{CHROMOSOME}-variantQC_filtered_FINAL.tsv.bgz", header=True)



