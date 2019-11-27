
'''
QC for 1 chromosome
author: Pavlos Antoniou
date: 22/07/19
'''

import os
import hail as hl
import sys
import json



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
    centromere_table = hl.import_bed("gs://interval-wgs/qc-files/Centromere_region_UCSC_GRCh38.bed", reference_genome='GRCh38', min_partitions = 250)

    #####################################################################
    ###################### INPUT DATA  ##############################
    #####################################################################
    # Give chromosome as input to program with chr prefix i.e chr1, chr2, chr3 etc
    CHROMOSOME = sys.argv[1]
    print(f"Reading {CHROMOSOME} mt")
    mt = hl.read_matrix_table(f"{BUCKET}/checkpoints/{CHROMOSOME}/{CHROMOSOME}_GT_fixed_NEW.mt")

    #mt = mt.key_rows_by('locus').distinct_by_row().key_rows_by('locus', 'alleles')

    print("Splitting mt and writing out split mt")
    mt_split = hl.split_multi(mt, keep_star=False, left_aligned=False)

    mt_split = mt_split.checkpoint(f"{BUCKET}/checkpoints/{CHROMOSOME}/{CHROMOSOME}-split-multi_checkpoint.mt", overwrite=True)
    print("Finished splitting and writing mt. ")

   