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
               "chr22",
               "chrX",
               "chrY"
               ]




if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    #Define the hail persistent storage directory
    hl.init(default_reference="GRCh38", tmp_dir=tmp_dir)
    fields_to_drop = ['PGT', 'PID']
    fields_to_drop_secondmt = ['ClippingRankSum', 'RAW_MQ']
    mt_chr1 = hl.read_matrix_table(f"{BUCKET}/checkpoints/chr1/chr1-split-multi_checkpoint.mt")
    mt_chr1=mt_chr1.drop(*fields_to_drop)
    info2 = mt_chr1.info.drop(*fields_to_drop_secondmt)
    mt_chr1 = mt_chr1.annotate_rows(info=info2)
    mt_chr1 = mt_chr1.annotate_entries(SB=mt_chr1.SB.map(lambda x: x))


    for CHROMOSOME in CHROMOSOMES:
        print(f"Reading chromosome {CHROMOSOME}")
        mt = hl.read_matrix_table(f"{BUCKET}/checkpoints/{CHROMOSOME}/{CHROMOSOME}-split-multi_checkpoint.mt")
        if CHROMOSOME!="chrX" and CHROMOSOME !="chrY":
            mt = mt.drop(*fields_to_drop)
            info2 = mt.info.drop(*fields_to_drop_secondmt)
            mt = mt.annotate_rows(info=info2)
            mt = mt.annotate_entries(SB=mt.SB.map(lambda x: x))
        mt_chr1 = mt_chr1.union_rows(mt)

    CHROMOSOME = "WGS_unfiltered"
    print("Finished merging successfully!!!!! ")


    #####################################################################
    ###################### QC ####################
    #####################################################################


    mt2 = hl.sample_qc(mt_chr1, name='sample_QC_Hail')
    mt3 = hl.variant_qc(mt2, name='variant_QC_Hail')
    print("Writing out")
    mt3 = mt3.checkpoint(
        f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-full-sampleqc-variantqc-UNFILTERED.mt", overwrite=True)
    mt3_cols = mt3.cols()
    mt3_cols.flatten().export(
        f"{BUCKET}/output-tables/{CHROMOSOME}/{CHROMOSOME}-sampleQC_UNFILTERED.tsv.bgz", header=True)

    mt3_rows = mt3.rows()
    mt3_rows.select(mt3_rows.variant_QC_Hail).flatten().export(
        f"{BUCKET}/output-tables/{CHROMOSOME}/{CHROMOSOME}-variantQC_UNFILTERED.tsv.bgz", header=True)


