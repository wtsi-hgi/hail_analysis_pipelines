'''
UKB gwas 
author: Pavlos Antoniou
date: 15/01/2020
'''

import os
import hail as hl
import pyspark
import json
import sys
from pathlib import Path
from bokeh.plotting import output_file, save




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
    tmp_dir = "hdfs://spark-master:9820"
    temp_dir = os.path.join(os.environ["HAIL_HOME"], "tmp")
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    #s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])

    phenotypes=hl.import_table("s3a://ukbb-plink/phenotypes-final.txt", delimiter='\t', types={'disease': hl.tfloat64})

    partitions=250
    #vcf=f"{temp_dir}/ukbb-plink/ukbb-plink-rename.vcf.bgz"
    #mt = hl.import_vcf(vcf, force_bgz=True,reference_genome='GRCh38') #, skip_invalid_loci=True)
    #if mt.n_partitions() > partitions:
    #    mt = mt.naive_coalesce(partitions)
    #mt= mt.checkpoint(f"{tmp_dir}/ukbb-plink/ukbb-plink.mt", overwrite=True)
    mt=hl.read_matrix_table(f"{tmp_dir}/ukbb-plink/ukbb-plink.mt")
    ja=phenotypes.key_by('samplename')
    mt = mt.annotate_cols(phenotype=ja[mt.s])
    mt= mt.checkpoint(f"{tmp_dir}/ukbb-plink/ukbb-plink_annotated.mt", overwrite=True)
    print("Run gwas")
    gwas = hl.linear_regression_rows(
        y=mt.phenotype.disease,
        x=mt.GT.n_alt_alleles(), covariates=[1.0], pass_through=[mt.rsid])
    gwas = gwas.checkpoint(f"{tmp_dir}/ukbb-plink/gwas/ukbb-plink-checkpoint", overwrite=True)


    print(f"Plotting ")
    p = hl.plot.manhattan(gwas.p_value, title="UKBB plink Manhattan plot")
    output_file(f"{tmp_dir}/ukbb-plink/gwas/plots/UKBB-manhattan.html")
    save(p)
        #hl.hadoop_copy(f"{i}.WGS-manhattan-{covariates[i]}.html", f"{BUCKET}/gwas/plots/")
    p = hl.plot.qq(gwas.p_value, title="UKBB plink QQ plot")
    output_file(f"{tmp_dir}/ukbb-plink/gwas/plots/UKBB-QQplot.html")
    save(p)
