import os
from datetime import datetime
import yaml
import hail as hl
import re
from bokeh.plotting import output_file, save
import json

project_root = os.path.dirname(os.path.dirname(__file__))
print(project_root)


snp_vqsr_threshold= -0.6647
indel_vqsr_threshold=-1.2537



BUCKET = "gs://interval-wgs"
#Define chromosome here
tmp_dir="/Users/pa10/Programming/google-code/google/tmp"


CHROMOSOMES = ["chr2",
               "chr3",
               "chr4",
               "chr5",
               "chr6",
               "chr7",
               "chr8",
               "chr9",
               "chr10",
               "chr11",
               "chr12",
               "chr13",
               "chr14",
               "chr15",
               "chr16",
               "chr17",
               "chr18",
               "chr19",
               "chr20",
               "chr21",
               "chr22"
               ]
               #"chrX",
               #"chrY"
               #]



covariates =[
    "neut_p_gwas_normalised",
    "eo_p_gwas_normalised",
    "mono_p_gwas_normalised",
   "lymph_p_gwas_normalised",
    "baso_p_gwas_normalised",
   "ret_p_gwas_normalised",
   "hlr_p_gwas_normalised",
   "hct_gwas_normalised",
   "pct_gwas_normalised",
   "hgb_gwas_normalised",
   "rbc_gwas_normalised",
   "wbc_gwas_normalised",
   "mpv_gwas_normalised",
   "plt_gwas_normalised",
   "rdw_gwas_normalised",
   "pdw_gwas_normalised",
   "mcv_gwas_normalised",
   "mch_gwas_normalised",
   "mchc_gwas_normalised",
   "ret_gwas_normalised",
   "hlr_gwas_normalised",
   "neut_gwas_normalised",
   "mono_gwas_normalised",
   "baso_gwas_normalised",
   "eo_gwas_normalised",
   "lymph_gwas_normalised",
   "irf_gwas_normalised",
   "myeloid_wbc_gwas_normalised",
   "gran_gwas_normalised",
   "eo_baso_sum_gwas_normalised",
   "neut_eo_sum_gwas_normalised",
   "baso_neut_sum_gwas_normalised",
   "gran_p_myeloid_wbc_gwas_normalised",
   "eo_p_gran_gwas_normalised",
   "neut_p_gran_gwas_normalised",
   "baso_p_gran_gwas_normalised"

]




if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    #Define the hail persistent storage directory
    hl.init(default_reference="GRCh38", tmp_dir=tmp_dir)


    mt_chr1 = hl.read_matrix_table(f"{BUCKET}/checkpoints/chr1/chr1-split-multi_checkpoint.mt")


    for CHROMOSOME in CHROMOSOMES:
        mt = hl.read_matrix_table(f"{BUCKET}/checkpoints/{CHROMOSOME}/{CHROMOSOME}-split-multi_checkpoint.mt")
        mt_chr1 = mt_chr1.union_rows(mt)

    CHROMOSOME = "WGS-autosomes"



    #####################################################################
    ###################### QC ####################
    #####################################################################


    mt2 = hl.sample_qc(mt_chr1, name='sample_QC_Hail')
    mt3 = hl.variant_qc(mt2, name='variant_QC_Hail')

    mt3 = mt3.checkpoint(
        f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-full-sampleqc-variantqc-UNFILTERED.mt", overwrite=True)
    mt3_cols = mt3.cols()
    mt3_cols.flatten().export(
        f"{BUCKET}/output-tables/{CHROMOSOME}/{CHROMOSOME}-sampleQC_UNFILTERED.tsv.bgz", header=True)

    mt3_rows = mt3.rows()
    mt3_rows.select(mt3_rows.variant_QC_Hail).flatten().export(
        f"{BUCKET}/output-tables/{CHROMOSOME}/{CHROMOSOME}-variantQC_UNFILTERED.tsv.bgz", header=True)


