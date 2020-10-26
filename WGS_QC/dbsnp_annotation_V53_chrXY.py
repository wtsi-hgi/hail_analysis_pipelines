import argparse
import sys
import logging
from pprint import pformat
import os
import argparse
import hail as hl
import pandas as pd
import numpy as np
import pyspark
import json
import sys
import re
from pathlib import Path
import logging
from typing import Any, Counter, List, Optional, Tuple, Union, Dict
import uuid
import json

os.environ['PYSPARK_PYTHON'] = sys.executable

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

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

tmp_dir = "hdfs://spark-master:9820/"
temp_dir = "file:///home/ubuntu/data/tmp"
plot_dir = "/home/ubuntu/data/tmp"
if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory

    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])
    n_partitions = 500


dbsnp_v153 = hl.read_matrix_table(
    f'{temp_dir}/intervalwgs/dbsnp_v53/dbsnp_v53_split.mt')


mt = hl.read_matrix_table(
    f'{temp_dir}/intervalwgs/WGS_final_march_2020_dbsnp_v53.mt')
mt_chrX = mt.filter_rows(mt.locus.contig == "chrX")
mt_chrX = hl.split_multi_hts(mt_chrX)
mt_chrX = mt_chrX.annotate_rows(rsid=dbsnp_v153.rows()[mt_chrX.row_key].rsid)
mt_chrX = mt_chrX.checkpoint(
    f"{tmp_dir}/chrX_dbsnp_v153_split.mt", overwrite=True)
mt_chrY = mt.filter_rows(mt.locus.contig == "chrY")
mt_chrY = hl.split_multi_hts(mt_chrY)
mt_chrY = mt_chrY.annotate_rows(rsid=dbsnp_v153.rows()[mt_chrY.row_key].rsid)
mt_chrY = mt_chrY.checkpoint(
    f"{tmp_dir}/chrY_dbsnp_v153_split.mt", overwrite=True)

hl.export_vcf(
    mt_chrX, f"{tmp_dir}/chrX_dbspn_v153.vcf.bgz")
hl.export_vcf(
    mt_chrY, f"{tmp_dir}/chrY_dbspn_v153.vcf.bgz")

mtX = mt_chrX.select_entries(mt_chrX.GT)
mtY = mt_chrY.select_entries(mt_chrY.GT)

hl.export_vcf(
    mtX, f"{tmp_dir}/chrX_dbspn_v153_GT_only.vcf.bgz")
hl.export_vcf(
    mtY, f"{tmp_dir}/chrY_dbspn_v153_GT_only.vcf.bgz")
