
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


def print_clusters(phenotype_group):
    print("**************************************************")
    nmr_dict={}
    for index,value in enumerate(phenotype_group):
        list1=value.collect()
        list2=[0 if v is None else 1 for v in list1]
        nmr_dict[nmr2[index]]=list2

    dict1={}
    namelist=[]
    for name, measurements in nmr_dict.items():

        tuple1=tuple(measurements)

        if tuple1 in dict1:
            # append the new number to the existing array at this slot
            dict1[tuple1].append(name)
            print("appended")
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
        
    nmr_dict={}
    all_groups=[nmr,somalogic_proteomics,olink_proteomics,fbc]
    for group in all_groups:
        print_clusters(group)
       
    print("Done")