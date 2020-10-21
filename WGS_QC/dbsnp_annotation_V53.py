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


dbsnp_vcf = "s3a://intervalwgs-qc/dbsnp_153.hg38.vcf.gz"
dbsnp_mt = hl.import_vcf(dbsnp_vcf, force_bgz=True, reference_genome='GRCh38',
                         skip_invalid_loci=True)
dbsnp_mt = dbsnp_mt.key_rows_by(
    'locus').distinct_by_row().key_rows_by('locus', 'alleles')
dbsnp_mt = hl.split_multi_hts(dbsnp_mt)
dbsnp_mt.write(f"{tmp_dir}/dbsnp_v53_split.mt", overwrite=True)
hl.export_vcf(dbsnp_mt, f"{tmp_dir}/dbsnp_v53_split.vcf.bgz", overwrite=True)
