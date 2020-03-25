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


def print_clusters(phenotype_group, phenotype_name,filename):
    print("**************************************************")
    nmr_dict={}
    #Create a dictionary with values a list of 1 if measurement in sample, 0 if not. 
    #key is the phenotype
    for index,value in enumerate(phenotype_group):
        list1=value.collect()
        list2=[0 if v is None else 1 for v in list1]
        nmr_dict[phenotype_name[index]]=list2

    dict1={}
    namelist=[]
    #convert each dictionary values in tuples
    #if this tuple is already in the dictionary (same samples)
    #append the phenotype name to the dictionary with key the tuple of measurements
    #if it's not create a new key-value pair with 
    for name, measurements in nmr_dict.items():

        tuple1=tuple(measurements)

        if tuple1 in dict1:
            # append the new number to the existing array at this slot
            dict1[tuple1].append(name)
        else:
            # create a new array in this slot
            dict1[tuple1] = [name]

    for key,value in dict1.items():
   # print(str(key) + '=>'+ str(value) + '=>' + str(len(value)))
        
        print( str(value))
        filename.write(str(value))
        filename.write("\n")

project="INT"
dataset="WGS"

nmr=[]
sysmex=[]
metabolon_metabolomics=[]
olink_inf=[]
olink_cvd2=[]
olink_cvd3=[]
olink_neu=[]
somalogic_proteomics=[]
fbc=[]
pcas=[]
nmr2=[]
sysmex2=[]
olink2_inf=[]
olink2_cvd2=[]
olink2_cvd3=[]
olink2_neu=[]
somalogic2=[]
fbc2=[]
pcas_names=[]
covariates_array=[]
covariates_names=[]
covariates_olinkf=[]
covariates_olinkf_names=[]
covariates_olinkcvd2=[]
covariates_olinkcvd2_names=[]
covariates_olinkcvd3=[]
covariates_olinkcvd3_names=[]
covariates_olinkneu=[]
covariates_olinkneu_names=[]

if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    #Define the hail persistent storage directory
    tmp_dir = "hdfs://spark-master:9820/"
    temp_dir = os.path.join(os.environ["HAIL_HOME"], "tmp")
    now= datetime.now()

    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38", log=temp_dir +f'/logfile-{now}.log')
    #s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])
     
    phenotypes = hl.import_table("s3a://intervalwgs-qc/phenotypes_wgs_wes_fbc_sysmex_nmr_olink_somalogic_50-wes-removed_05-03-2020.csv", impute=True, delimiter=',')
    #dbsnp = hl.import_vcf("s3a://intervalwgs-qc/GCF_000001405.38.fixed.vcf.gz", force_bgz=True, reference_genome='GRCh38',skip_invalid_loci=True)
    #dbsnp_rows = dbsnp.rows()

    CHROMOSOME="WGS"

    f = open(f"{temp_dir}/scripts/hail-pipelines-internal/gwas/phenotype_lists/olink_phenotype_names.txt", "w")
    mt = hl.read_matrix_table(f"{temp_dir}/intervalwgs/WGS_final_february_2020_updated_rsID.mt")
    
    print("Number of initial variants:")
    print(mt.count())
    
    print("Annotating matrixtable with phenotypes:")
    phenotypes=phenotypes.key_by('ID')
    mt = mt.annotate_cols(phenotype=phenotypes[mt.s])
    

    print("Grouping the phenotypes into lists:")
    ph1=list(mt.phenotype)
    for pheno in ph1:
        if pheno =='BWA' or pheno =='Study':
            covariates_array.append(hl.float64(mt.phenotype[pheno]))
            covariates_names.append(pheno)
        if pheno.startswith('somalogic'):
            somalogic_proteomics.append(mt.phenotype[pheno])
            somalogic2.append(pheno)
        elif pheno.startswith('nmr'):
            nmr.append(mt.phenotype[pheno])
            nmr2.append(pheno)
        elif pheno.startswith('sysmex'):
            sysmex.append(mt.phenotype[pheno])
            sysmex2.append(pheno)
        elif pheno.startswith('olinkinf'):
            if pheno == 'olinkf_plate':
                covariates_olinkf.append(hl.float64(mt.phenotype[pheno]))
                covariates_olinkf_names.append(pheno)
            else:
                olink_inf.append(mt.phenotype[pheno])
                olink2_inf.append(pheno)
        elif pheno.startswith('olinkcvd2'):
            if pheno == 'olinkcvd2_plate':
                covariates_olinkcvd2.append(hl.float64(mt.phenotype[pheno]))
                covariates_olinkcvd2_names.append(pheno)
            else:
                olink_cvd2.append(mt.phenotype[pheno])
                olink2_cvd2.append(pheno)
        elif pheno.startswith('olinkcvd3'):
            if pheno == 'olinkcvd3_plate':
                covariates_olinkcvd3.append(hl.float64(mt.phenotype[pheno]))
                covariates_olinkcvd3_names.append(pheno)
            else:
                olink_cvd3.append(mt.phenotype[pheno])
                olink2_cvd3.append(pheno)
        elif pheno.startswith('olinkneu'):
            if pheno == 'olinkneu_plate':
                covariates_olinkneu.append(hl.float64(mt.phenotype[pheno]))
                covariates_olinkneu_names.append(pheno)
            else:
                olink_neu.append(mt.phenotype[pheno])
                olink2_neu.append(pheno)
        elif pheno.startswith('fbc'):
            fbc.append(mt.phenotype[pheno])
            fbc2.append(pheno)
        elif pheno.startswith('PC'):
            pcas.append(mt.phenotype[pheno])
            pcas_names.append(pheno)

    covariates_array=pcas+covariates_array
    covariates_names=pcas_names+covariates_names

    #all_groups=[nmr,somalogic_proteomics,sysmex,olink_inf, olink_cvd2,olink_cvd3,olink_neu, fbc]
    #all_names=[nmr2,somalogic2,somalogic2, sysmex2, olink2_inf, olink2_cvd2,olink2_cvd3,olink2_neu, fbc2]

    all_groups=[olink_inf, olink_cvd2,olink_cvd3,olink_neu]
    all_names=[olink2_inf, olink2_cvd2,olink2_cvd3,olink2_neu]
    for list_name in all_names:
        for phenotype in list_name:
            print(phenotype)
            f.write(phenotype)
            f.write("\n")
    print("Done")