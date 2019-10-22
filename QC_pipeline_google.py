
'''
QC for 1 chromosome
author: Pavlos Antoniou
date: 22/07/19
'''

import os
import hail as hl
import sys
import json

CHROMOSOME="chr20"

project_root = os.path.dirname(os.path.dirname(__file__))
print(project_root)



snp_vqsr_threshold= -0.6647
indel_vqsr_threshold=-1.2537



BUCKET = "gs://interval-wgs"
#Define chromosome here
tmp_dir="/Users/pa10/Programming/google-code/google/tmp"

if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    #Define the hail persistent storage directory
    hl.init(default_reference="GRCh38", tmp_dir=tmp_dir)


    #1. Import required files
    print("1. import required tables in hail for the project.")
    VQSLOD_snps = hl.import_table("gs://interval-wgs/qc-files/VQSLOD_snps.bgz",
                                  types={"Locus": "locus<GRCh38>", "VQSLOD": hl.tfloat64})
    VQSLOD_indels = hl.import_table("gs://interval-wgs/qc-files/VQSLOD_indels.bgz",
                                    types={"Locus": "locus<GRCh38>", "VQSLOD": hl.tfloat64})
    sample_QC_nonHail = hl.import_table("gs://interval-wgs/qc-files/INTERVAL_WGS_Sample_QC_04-09-2019.txt", impute=True)

    #####################################################################
    ###################### INPUT DATA  ##############################
    #####################################################################
    # Give chromosome as input to program with chr prefix i.e chr1, chr2, chr3 etc
    #CHROMOSOME = sys.argv[1]
    print(f"Reading {CHROMOSOME} mt")
    mt = hl.read_matrix_table(f"{BUCKET}/{CHROMOSOME}.mt")

    print("Splitting mt and writing out split mt")
    mt_split = hl.split_multi_hts(mt, keep_star=False)

    mt_split = mt_split.checkpoint(f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-split-multi.mt", overwrite=True)
    print("Finished splitting and writing mt. ")

    #####################################################################
    ###################### UNFILTERED SAMPLE AND VARIANT QC #############
    #####################################################################

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
    panda_df_unfiltered_table.export(f"{BUCKET}/output-tables/{CHROMOSOME}/{CHROMOSOME}_sampleQC_unfiltered.tsv.bgz", header=True)

    # Variant QC
    mt2 = hl.variant_qc(mt_sqc2_unfiltered, name='variant_QC_Hail')

    print('Exporting variant qc pandas table to disk')
    mt_rows = mt2.rows()
    mt_rows.select(mt_rows.variant_QC_Hail).flatten().export(
        f"{BUCKET}/output-tables/{CHROMOSOME}/{CHROMOSOME}_variantQC_unfiltered.tsv.bgz",
        header=True)

    #####################################################################
    ###################### START FILTERING ##############################
    #####################################################################

    ######## VQSR filtering
    print("Annotation VQSLOD snp and indel scores")
    mt = mt.annotate_rows(VQSLOD_SNP=VQSLOD_snps.key_by("Locus")[mt.locus])
    mt = mt.annotate_rows(VQSLOD_INDEL=VQSLOD_indels.key_by("Locus")[mt.locus])

    print("Filtering on VQSLOD scores")
    vqslod_filtered = (
        hl.case()
            .when((mt.Variant_Type == "SNP"), (mt.VQSLOD_SNP.VQSLOD >= snp_vqsr_threshold))
            .when((mt.Variant_Type == "INDEL"),
                  (mt.VQSLOD_INDEL.VQSLOD >= indel_vqsr_threshold))
            .default(False)  # remove everything else
    )

    mt_vqslod_filtered = mt.filter_rows(vqslod_filtered)
    print("Finished filtering. Writing out matrixtable...")

    mt1 = mt_vqslod_filtered.checkpoint(f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}_vqslod_filtered_checkpoint.mt",
                                        overwrite=True)
    print("Finished writing vqslod filtered matrixtable")

    ########### Sample QC filtering
    print("Filtering on sample qc")
    mt_sqc1_filtered = mt1.filter_cols(
        (mt1.sample_QC_nonHail.PASS_Depth == 1) &
        (mt1.sample_QC_nonHail.PASS_ID == 1) &
        (mt1.sample_QC_nonHail.PASS_Median_FreeMix == 1) &
        (mt1.sample_QC_nonHail.PASS_NRD == 1) &
        (mt1.sample_QC_nonHail.PASS_SampleSwap == 1) &
        (mt1.sample_QC_nonHail.PASS_Sex == 1) &
        (mt1.sample_QC_nonHail.PASS_DUP == 1)
    )
    print("Writing out filtered sample qc checkpoint")
    mt_sqc2 = hl.sample_qc(mt_sqc1_filtered, name='sample_QC_Hail')

    mt_sqc2 = mt_sqc2.checkpoint(f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}_sampleQC_filtered_checkpoint.mt",
                                 overwrite=True)
    filter_sampleqc_table = mt_sqc2.cols().flatten()
    filter_sampleqc_table.export(f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-sampleQC_filtered.tsv.bgz", header=True)

    ##### ADD allele bias script here

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

    mt_sqc2_GT = mt_sqc2.filter_entries(filter_condition_ab, keep=False)

    ## Saving on s3
    # mt_sqc2 = mt_sqc2_GT.checkpoint('s3a://interval-wgs/checkpoints_new/chr20_sampleQC_step1_filtered_allele_balance_checkpoint.mt', _read_if_exists = True)

    ## Saving on local volume
    mt_sqc2 = mt_sqc2_GT.checkpoint(
        f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}_sampleQC_step1_filtered_allele_balance_checkpoint.mt",
        overwrite=True)

    ######################
    ##### VARIANT qc
    ######################
    print("Variant qc:")
    mt_sqc_vqc = hl.variant_qc(mt_sqc2, name='variant_QC_Hail')
    mt_sqc_vqc_filtered = mt_sqc_vqc.filter_rows(
        (mt_sqc_vqc.variant_QC_Hail.call_rate >= 0.98) &
        (mt_sqc_vqc.variant_QC_Hail.p_value_hwe >= 10 ** -6))

    #####################################################################
    ###################### FINAL QC AFTER FILTERING  ####################
    #####################################################################

    fields_to_drop = ['variant_QC_Hail', 'sample_QC_Hail']

    mt1 = mt_sqc_vqc_filtered.drop(*fields_to_drop)

    mt2 = hl.sample_qc(mt1, name='sample_QC_Hail')
    mt3 = hl.variant_qc(mt2, name='variant_QC_Hail')

    mt3 = mt_sqc_vqc_filtered.checkpoint(
        f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-full-sampleqc-variantqc-filtered-FINAL.mt", overwrite=True)
    mt3_rows = mt3.rows()
    mt3_rows.select(mt3_rows.variant_QC_Hail).flatten().export(
        f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-variantQC_sampleQC_filtered_FINAL.tsv.bgz", header=True)

