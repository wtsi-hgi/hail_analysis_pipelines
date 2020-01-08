
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
nmr_samples=11410
olink_inf_samples=2984
olink_cvd2_samples=3015
olink_cvd3_samples=3044
olink_neu_samples=2859
somalogic_samples=55
fbc_samples=11096
sysmex_samples=11754

def print_clusters(phenotype_group, phenotype_name, group_samples):
    print("**************************************************")
    nmr_dict={}
    #Create a dictionary with values a list of 1 if measurement in sample, 0 if not. 
    #key is the phenotype
    i=1
    print("id\tPhenotype\tsamples\tTotalsamples\tFraction\tMissingess\tSamples_missing")
    for phenotype,name in zip(phenotype_group,phenotype_name):
        numofsamples=mt.aggregate_cols(hl.agg.count_where(hl.is_defined(phenotype)))
        mt2=mt.filter_cols(hl.is_defined(phenotype),keep=False)
        samples=mt2.s.collect()
        #print(samples)
        fraction=(int(numofsamples)*100/group_samples)
        missingness=100-fraction
        print(f'{i}\t{name}\t{numofsamples}\t{group_samples}\t{fraction:.3f}\t{missingness:.3f}\t{samples}')
        i+=1
    
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
        elif pheno.startswith('sysmex'):
            sysmex.append(mt.phenotype[pheno])
            sysmex2.append(pheno)
        elif pheno.startswith('olinkinf'):
            olink_inf.append(mt.phenotype[pheno])
            olink2_inf.append(pheno)
        elif pheno.startswith('olinkcvd2'):
            olink_cvd2.append(mt.phenotype[pheno])
            olink2_cvd2.append(pheno)
        elif pheno.startswith('olinkcvd3'):
            olink_cvd3.append(mt.phenotype[pheno])
            olink2_cvd3.append(pheno)
        elif pheno.startswith('olinkneu'):
            olink_neu.append(mt.phenotype[pheno])
            olink2_neu.append(pheno)
        elif pheno.startswith('fbc'):
            fbc.append(mt.phenotype[pheno])
            fbc2.append(pheno)
        elif pheno.startswith('PC'):
            pcas.append(mt.phenotype[pheno])
            
        
        

    #all_groups=[nmr,somalogic_proteomics,olink_proteomics,fbc]
    #all_names=[nmr2,somalogic2,olink2,fbc2]

    all_groups=[nmr,sysmex,
    olink_inf,olink_cvd2,olink_cvd3,olink_neu,fbc]

    all_names=[nmr2,sysmex2,
    olink2_inf,olink2_cvd2,olink2_cvd3,olink2_neu,fbc2]

    all_samples=[nmr_samples,sysmex_samples,
   olink_inf_samples,olink_cvd2_samples,olink_cvd3_samples,
   olink_neu_samples,fbc_samples]
    
    for group,name,samples in zip(all_groups,all_names,all_samples):
        print_clusters(group,name,samples)
       
    print("Done")
