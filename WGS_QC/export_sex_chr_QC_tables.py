
'''
QC for 1 chromosome
author: Pavlos Antoniou
date: 10/01/2020
'''

import os
import hail as hl
import pyspark
import json
import sys
from pathlib import Path


project_root = Path(__file__).parent.parent.parent
print(project_root)

s3credentials = os.path.join(project_root, "config_files/s3_credentials.json")
print(s3credentials)

storage = os.path.join(project_root, "config_files/storage.json")

thresholds = os.path.join(project_root, "config_files/thresholds.json")

with open(f"{s3credentials}", 'r') as f:
    credentials = json.load(f)

with open(f"{storage}", 'r') as f:
    storage = json.load(f)

with open(f"{thresholds}", 'r') as f:
    thresholds = json.load(f)


if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory
    tmp_dir = "hdfs://spark-master:9820/"
    temp_dir = os.path.join(os.environ["HAIL_HOME"], "tmp")
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])

    # 1. Import required files
    print("1. import required tables in hail for the project.")
    VQSLOD_snps = hl.import_table(storage["intervalwgs"]["s3"]["vqsr_snp"],
                                  types={"Locus": "locus<GRCh38>", "VQSLOD": hl.tfloat64})
    VQSLOD_indels = hl.import_table(storage["intervalwgs"]["s3"]["vqsr_indel"],
                                    types={"Locus": "locus<GRCh38>", "VQSLOD": hl.tfloat64})
    sample_QC_nonHail = hl.import_table(
        storage["intervalwgs"]["s3"]["sampleQC_non_hail"], impute=True)

    centromere_table = hl.import_bed(
        storage["intervalwgs"]["s3"]["centromere"], reference_genome='GRCh38', min_partitions=250)

    #####################################################################
    ###################### INPUT DATA  ##############################
    #####################################################################
    # Give chromosome as input to program with chr prefix i.e chr1, chr2, chr3 etc
    mtX = hl.read_matrix_table(
        f"{temp_dir}/intervalwgs/sex_chromosomes/chrX_final_AC_corrected_filtered_variant_QC_ordered_FINAL.mt")
    mtY = hl.read_matrix_table(
        f"{temp_dir}/intervalwgs/sex_chromosomes/chrY_final_AC_corrected_filtered_variant_QC_ordered_FINAL.mt")

    mtX_cols = mtX.cols()
    mtX_cols.flatten().export(
        f"{tmp_dir}/intervalwgs/sex_chromosomes/chrX-sampleQC_filtered_FINAL.tsv.bgz", header=True)

    mtX_rows = mtX.rows()
    mtX_rows.select(mtX_rows.variant_QC_Hail).flatten().export(
        f"{tmp_dir}/intervalwgs/sex_chromosomes/chrX-variantQC_filtered_FINAL.tsv.bgz", header=True)

    mtY_cols = mtY.cols()
    mtY_cols.flatten().export(
        f"{tmp_dir}/intervalwgs/sex_chromosomes/chrY-sampleQC_filtered_FINAL.tsv.bgz", header=True)

    mtY_rows = mtY.rows()
    mtY_rows.select(mtY_rows.variant_QC_Hail).flatten().export(
        f"{tmp_dir}/intervalwgs/sex_chromosomes/chrY-variantQC_filtered_FINAL.tsv.bgz", header=True)
