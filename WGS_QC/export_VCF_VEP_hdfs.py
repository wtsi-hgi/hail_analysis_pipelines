
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
    CHROMOSOME = "chrX"
    mt_wgs = hl.read_matrix_table(
        f"{temp_dir}/intervalwgs/WGS_final_march_2020_dbsnp_v53.mt")
    mt = mt_wgs.filter_rows(mt_wgs.locus.contig == CHROMOSOME)
    mt1 = mt.select_entries()
    mt_fin = mt1.filter_cols(mt1['s'] == 'sample')
    hl.export_vcf(
        mt_fin, f"{tmp_dir}/VCFs/{CHROMOSOME}/{CHROMOSOME}_for_VEP.vcf.bgz")
    CHROMOSOME = "chrY"
    mt_wgs = hl.read_matrix_table(
        f"{temp_dir}/intervalwgs/WGS_final_march_2020_dbsnp_v53.mt")
    mt = mt_wgs.filter_rows(mt_wgs.locus.contig == CHROMOSOME)
    mt1 = mt.select_entries()
    mt_fin = mt1.filter_cols(mt1['s'] == 'sample')
    hl.export_vcf(
        mt_fin, f"{tmp_dir}/VCFs/{CHROMOSOME}/{CHROMOSOME}_for_VEP.vcf.bgz")
