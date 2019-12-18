import os
from datetime import datetime
import yaml
import hail as hl
import re
from bokeh.plotting import output_file, save
import json

project_root = os.path.dirname(os.path.dirname(__file__))
print(project_root)


snp_vqsr_threshold= -0.6647
indel_vqsr_threshold=-1.2537



BUCKET = "gs://interval-wgs"
#Define chromosome here
tmp_dir="/Users/pa10/Programming/google-code/google/tmp"

""" 
print(len(nmr))
print(len(metabolon_metabolomics))
print(len(olink_proteomics))
print(len(somalogic_proteomics))
print(len(fbc))
print(len(pcas))
print(len(nmr)+len(metabolon_metabolomics)+len(olink_proteomics)+len(somalogic_proteomics)+len(fbc)+len(pcas))
print(len(ph1))
230
0
372
4035
36
10
4683
4683 
"""


nmr=[]
metabolon_metabolomics=[]
olink_proteomics=[]
somalogic_proteomics=[]
fbc=[]
pcas=[]



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
        elif pheno.startswith('olink'):
            olink_proteomics.append(mt.phenotype[pheno])
        elif pheno.startswith('fbc'):
            fbc.append(mt.phenotype[pheno])
        elif pheno.startswith('PC'):
            pcas.append(mt.phenotype[pheno])
        
    

    print("Linear regression")
    #y=[nmr1[0],fbc1[0], op[0],sp[0]],
    #y=[nmr1[0],fbc1[0], op[0],sp[0]],
    # TIM NOTE: I changed the below to show how linear_regression_rows can process groups of phenotypes in a vectorized way
    gwas = hl.linear_regression_rows(
        y=[nmr],
        x=mt.GT.n_alt_alleles(), covariates=[1.0]+pcas, pass_through=[mt.rsid])

    # gwas = hl.linear_regression_rows(y=[[mt_filtered.sample_qc_and_phenotype.wbc_gwas_normalised, ...], [family2...]],
    print("Linear regression CHECKPOINT")
    # TIM NOTE: checkpoint here to prevent multiple execution (write to a file, read that file)
    gwas = gwas.checkpoint(f"{BUCKET}/gwas/{CHROMOSOME}-gwasfbc-checkpoint-nmr", overwrite=True)
    gwas=gwas.annotate(nmr_phenotypes=nmr)
    gwas.export(f"{BUCKET}/gwas/gwas-{CHROMOSOME}-export-somalogic_p_value_0.05.tsv.bgz", header=True)

    # gwas = gwas.checkpoint(s3"{tmp_dir}/gwas_wbc_chr19_checkpoint.mt")
    print("Number of variants in gwas table:")
    print("gwas.count()")
    #gwas1=gwas.filter(gwas.p_value[0].any(lambda x: x < 0.05), keep=True)
    
    #gwas1 = gwas1.checkpoint(f"{BUCKET}/gwas/{CHROMOSOME}-gwasfbc-checkpoint-", overwrite=True)
    #print("Number of varians after filtering by keeing only variants that have any phenotype with p-value less than < 0.05")
    #print("Number of varians after filtering by keeing only variants that have any phenotype with p-value less than < 0.04")
    #print(gwas2.count())
    ##print(gwas1.count())
    #print("Number of varians after filtering by keeing only variants that have any phenotype with p-value less than < 0.03")
    #print(gwas3.count())
    #print("Linear regression output table")
    #gwas1.export(f"{BUCKET}/gwas/gwas-{CHROMOSOME}-export-somalogic_p_value_0.05.tsv.bgz", header=True)

    print("Plotting")

    #for i in range(len(nmr)):
       # print(f"Plotting {i}:{nmr[i]}")
       # p = hl.plot.manhattan(gwas.p_value[i], title=f"Interval WGS GWAS Manhattan Plot: {nmr[i]}")
      #  output_file(f"{i}.WGS-manhattan-{nmr[i]}.html")
       # save(p)
       # hl.hadoop_copy(f"{i}.WGS-manhattan-{nmr[i]}.html", f"{BUCKET}/gwas/plots/")

       # p = hl.plot.qq(gwas.p_value[i], title=f"Interval WGS GWAS QQ Plot: {covariates[i]}")
        # export_png(p,f"{BUCKET}/output-tables/wgs-qq-{covariates[i]}.html")
       # output_file(f"{i}.WGS-qq-{covariates[i]}.html")
       # save(p)
       # hl.hadoop_copy(f"{i}.WGS-qq-{covariates[i]}.html", f"{BUCKET}/gwas/plots/")

