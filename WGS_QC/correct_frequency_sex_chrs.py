'''
Coorect AC and AN for chrX and chrY
author: Pavlos Antoniou
date: 10/01/2020
'''

import os
import hail as hl
import pyspark
import json
import sys
from pathlib import Path



project_root=Path(__file__).parent.parent.parent
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
    tmp_dir = "hdfs://spark-master:9820/"
    temp_dir = os.path.join(os.environ["HAIL_HOME"], "tmp")
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    #s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])

    mt_wgs=hl.read_matrix_table(f"{temp_dir}/intervalwgs/WGS_after_surgery-full-sampleqc-variantqc-FILTERED.mt")
    sample_QC_nonHail = hl.import_table(storage["intervalwgs"]["s3"]["sampleQC_non_hail"], impute=True)
    mt=mt_wgs.annotate_cols(sex=sample_QC_nonHail.key_by("ID")[mt_wgs.s].SEX)
    mt=mt.annotate_rows(call_stats= hl.agg.call_stats(mt.GT, mt.alleles))
    mt=mt.annotate_rows(final_AC= mt.call_stats.AC)
    mt=mt.annotate_rows(final_AN= mt.call_stats.AN)
    mt=mt.annotate_rows(frequency= hl.float64(mt.final_AC[1]/mt.final_AN))
    #Select chrX only
    mtX=mt.filter_rows(mt.locus.contig == "chrX",keep=True)
    mt_rest=mt.filter_rows(mt.locus.contig == "chrX",keep=False)
    mtX=annotate_allele_counts_chrX(mtX)
    mt=mt_rest.union_rows(mtX)
    mtY=mt.filter_rows(mt.locus.contig == "chrY",keep=True)
    mt_rest=mt.filter_rows(mt.locus.contig == "chrY",keep=False)
    mtY=annotate_allele_counts_chrY(mtY)
    mt=mt_rest.union_rows(mtY)
    mt=mt.checkpoint(f"{temp_dir}/intervalwgs/WGS_final_january_2020.mt", overwrite=True)