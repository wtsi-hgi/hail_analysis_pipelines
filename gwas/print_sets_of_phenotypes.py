
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
somalogic2=[]
olink2=[]
fbc2=[]

def print_clusters(phenotype_group, phenotype_name):
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
        print( str(value) + '=>' + str(len(value)))

    print("**************************************************")


if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    #Define the hail persistent storage directory
    hl.init(default_reference="GRCh38", tmp_dir=tmp_dir)

    phenotypes = hl.import_table("gs://interval-wgs/gwas/phenotypes_pc_fbc_nmr_olink_somalogic_16-12-1019.csv", impute=True, delimiter=',')

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
            somalogic2.append(pheno)
        elif pheno.startswith('nmr'):
            nmr.append(mt.phenotype[pheno])
            nmr2.append(pheno)
        elif pheno.startswith('olink'):
            olink_proteomics.append(mt.phenotype[pheno])
            olink2.append(pheno)
        elif pheno.startswith('fbc'):
            fbc.append(mt.phenotype[pheno])
            fbc2.append(pheno)
        elif pheno.startswith('PC'):
            pcas.append(mt.phenotype[pheno])
        

    #all_groups=[nmr,somalogic_proteomics,olink_proteomics,fbc]
    #all_names=[nmr2,somalogic2,olink2,fbc2]

    all_groups=[olink_proteomics,fbc]
    all_names=[olink2,fbc2]
    for group,name in zip(all_groups,all_names):
        print_clusters(group,name)
       
    print("Done")