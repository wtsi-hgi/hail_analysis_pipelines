
'''
gwas locally 
author: Pavlos Antoniou
date: 10/01/2020
'''

import os
import hail as hl
import pyspark
import json
import sys
from pathlib import Path
from bokeh.plotting import output_file, save, show
from datetime import datetime
import ast
from hail.utils.java import Env


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


project = "INT"
dataset = "WGS"

nmr = []
sysmex = []
metabolon_metabolomics = []
olink_inf = []
olink_cvd2 = []
olink_cvd3 = []
olink_neu = []
somalogic_proteomics = []
fbc = []
pcas = []
nmr2 = []
sysmex2 = []
olink2_inf = []
olink2_cvd2 = []
olink2_cvd3 = []
olink2_neu = []
somalogic2 = []
fbc2 = []
pcas_names = []
covariates_array = []
covariates_names = []


def GWAS_for_nmr(mt):
    pcas = []
    pcas_names = []
    covariates_array = []
    covariates_names = []
    ph1 = list(mt.phenotype)

    for pheno in ph1:
        if pheno == 'BWA' or pheno == 'Study' or pheno == "Recruitment_Centre":
            covariates_array.append(hl.float64(mt.phenotype[pheno]))
            covariates_names.append(pheno)
        if pheno.startswith('PC'):
            pcas.append(mt.phenotype[pheno])
            pcas_names.append(pheno)

    covariates_array = pcas+covariates_array
    covariates_names = pcas_names+covariates_names
    return covariates_array


def GWAS_for_olink(mt, groupname):
    pcas = []
    pcas_names = []
    covariates_array = []
    covariates_names = []
    plate = ""
    ph1 = list(mt.phenotype)
    if groupname == 'olinkinf':
        plate = 'olinkinf_plate'
    elif groupname == 'olinkcvd2':
        plate = 'olinkcvd2_plate'
    elif groupname == 'olinkcvd3':
        plate = 'olinkcvd3_plate'
    elif groupname == 'olinkneu':
        plate = 'olinkneu_plate'

    for pheno in ph1:

        if pheno == 'Study_BWA' or pheno == "Recruitment_Centre" or pheno == plate:
            covariates_array.append(hl.float64(mt.phenotype[pheno]))
            covariates_names.append(pheno)

        if pheno.startswith('PC'):
            pcas.append(mt.phenotype[pheno])
            pcas_names.append(pheno)
    covariates_array = pcas+covariates_array
    covariates_names = pcas_names+covariates_names

    print("The covariates are:")
    print(covariates_names)
    return covariates_array


if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory
    tmp_dir = "hdfs://spark-master:9820/"
    temp_dir = "file:///home/ubuntu/data"
    now = datetime.now()

    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])
    hadoop_config = Env().hc().sc._jsc.hadoopConfiguration()
    hadoop_config.set("fs.s3a.connection.maximum", "1000")
    gwas = hl.read_table(
        f'{tmp_dir}/gwas/INT-WGS-gwas-olinkcvd2-olinkcvd2_bmp6___p22004.tsv.bgz')
    index = 1
    pheno_name = "olinkcvd2_bmp6___p22004"
    print(f"Plotting manhattan {index}:{pheno_name}")
    p = hl.plot.manhattan(gwas.p_value, title=f"{pheno_name} GWAS")
    output_file(
        f"{temp_dir}/gwas/{project}-{dataset}-{pheno_name}-manhattan.html", mode='inline')
    save(p)
    print(f"Plotting QQ plot for {index} - {pheno_name}")
    q = hl.plot.qq(gwas.p_value, collect_all=False,
                   n_divisions=100, title=f"{pheno_name} QQ plot")
    output_file(
        f"{temp_dir}/gwas/{project}-{dataset}-{pheno_name}-QQplot.html", mode='inline')
    save(q)
