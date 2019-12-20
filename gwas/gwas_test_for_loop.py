
'''
author: Pavlos Antoniou
date: 22/07/19
'''

import os
import hail as hl
import pyspark
import json
import sys
from pathlib import Path


BUCKET = "gs://interval-wgs"
#Define chromosome here
tmp_dir="/Users/pa10/Programming/google-code/google/tmp"



nmr=[]
metabolon_metabolomics=[]
olink_proteomics=[]
somalogic_proteomics=[]
fbc=[]
pcas=[]
nmr2=[]

if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    #Define the hail persistent storage directory
    tmp_dir = os.path.join(os.environ["HAIL_HOME"], "tmp")
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    #s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])
    
    phenotypes = hl.import_table("gs://interval-wgs/gwas/phenotypes_pc_fbc_nmr_olink_somalogic_09-12-1019.csv", impute=True, delimiter=',')

    sample_QC_nonHail = hl.import_table("gs://interval-wgs/qc-files/INTERVAL_WGS_Sample_QC_04-09-2019.txt", impute=True)
    gws_gwa_map = hl.import_table(f"{BUCKET}/qc-files/WGS-2-GWA_omicsMap.txt", impute=True)


    print('Joining annotations')
    ja=phenotypes.key_by('ID')

    

    #after having merged chromosomes and done new sample and variant qc with merge_matrixtables_FINAL.py

    CHROMOSOME="WGS-autosomes"
    mt = hl.read_matrix_table(f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-full-sampleqc-variantqc-filtered-FINAL.mt")
    
    print("Number of initial variants:")
    print(mt.count())
    
    print("Annotating matrixtable with fbc:")
    mt = mt.annotate_cols(phenotype=ja[mt.s])

    print("Grouping the phenotypes into lists:")
    ph1=list(mt.phenotype)
    for pheno in ph1:
        if pheno.startswith('somalogic'):
            somalogic_proteomics.append(mt.phenotype[pheno])
        elif pheno.startswith('nmr'):
            nmr.append(mt.phenotype[pheno])
            nmr2.append(pheno)
        elif pheno.startswith('olink'):
            olink_proteomics.append(mt.phenotype[pheno])
        elif pheno.startswith('fbc'):
            fbc.append(mt.phenotype[pheno])
        elif pheno.startswith('PC'):
            pcas.append(mt.phenotype[pheno])
        
    

    print("Linear regression")
    for index,value in enumerate(nmr[0:10]):
        print(nmr2[index])
        gwas = hl.linear_regression_rows(
            y=value,
            x=mt.GT.n_alt_alleles(), covariates=[1.0]+pcas[0:3], pass_through=[mt.rsid])
        gwas = gwas.checkpoint(f"{BUCKET}/gwas/gwas{i}-test.table", overwrite=True)
        gwas.export(f"{BUCKET}/gwas/gwas-{i}_test_loop.tsv.bgz", header=True)

    print("Done")