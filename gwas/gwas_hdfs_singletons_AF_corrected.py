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

if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    #Define the hail persistent storage directory
    tmp_dir = "hdfs://spark-master:9820/"
    temp_dir = os.path.join(os.environ["HAIL_HOME"], "tmp")
    now= datetime.now()

    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    #s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])
     
    phenotypes = hl.import_table("s3a://intervalwgs-qc/phenotypes_wgs_wes_fbc_sysmex_nmr_olink_somalogic_50-wes-removed_05-03-2020.csv", impute=True, delimiter=',')
    #dbsnp = hl.import_vcf("s3a://intervalwgs-qc/GCF_000001405.38.fixed.vcf.gz", force_bgz=True, reference_genome='GRCh38',skip_invalid_loci=True)
    #dbsnp_rows = dbsnp.rows()

    CHROMOSOME="WGS"
    mt = hl.read_matrix_table(f"{temp_dir}/intervalwgs/WGS_final_march_2020_dbsnp_v53.mt")
    
    print("Number of initial variants:")
    print(mt.count())
    
    print("Annotating matrixtable with phenotypes:")
    phenotypes=phenotypes.key_by('ID')
    mt = mt.annotate_cols(phenotype=phenotypes[mt.s])


    #########################
    working_pheno_group="nmr"
    ########################
    print("Grouping the phenotypes into lists:")
    ph1=list(mt.phenotype)
    for pheno in ph1:
        if working_pheno_group == "nmr":
            if pheno =='BWA' or pheno =='Study' or pheno == "Recruitment_Centre":
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
            pcas_names.append(pheno)

    covariates_array=pcas+covariates_array
    covariates_names=pcas_names+covariates_names
    #make sure the input file has no empty line at the end 

    
        

    with open(f"{temp_dir}/scripts/hail-pipelines-internal/hail_analysis_pipelines/gwas/nmr_phenotypes.txt", 'r') as f:
        phenotypes_to_run=[ast.literal_eval(line.strip()) for line in f]

    nmr_new=[]
    nmr2_new=[]
    print("Linear regression")

    for group in phenotypes_to_run:

        print("Original group")
        print(group)
        for pheno in ph1:
            if pheno.startswith('nmr'):
                if pheno in group:
                    nmr_new.append(mt.phenotype[pheno])
                    nmr2_new.append(pheno)
        print("No running gwas with these phenotypes:")
        print(nmr2_new)
        mt = mt.annotate_rows(pheno_call_stats = hl.agg.filter(hl.is_defined(mt.phenotype[pheno]), hl.agg.call_stats(mt.GT, mt.alleles)))
        mt = mt.annotate_rows(pheno_n_het = hl.agg.filter(hl.is_defined(mt.phenotype[pheno]), hl.agg.count_where(mt.GT.is_het())))
        ## Applying MAF filter -- basically removing singletones and alternative AC = 0
        mt = mt.filter_rows((mt.pheno_call_stats.AC[0] != 1) &
                            (mt.pheno_call_stats.AC[1] >= 2)
                           )

        gwas = hl.linear_regression_rows(
            y=nmr_new[0],
            x=mt.GT.n_alt_alleles(), covariates=[1.0]+covariates_array, pass_through=[mt.rsid])
        fields_to_drop = ['sum_x', 'y_transpose_x','t_stat' ]
        gwas_table=gwas.drop(*fields_to_drop)
        gwas_table=gwas_table.annotate(nmr_phenotypes=nmr2_new[0])
        gwas_table=gwas_table.annotate(REF=gwas_table.alleles[0])
        gwas_table=gwas_table.annotate(ALT=gwas_table.alleles[1])
        #gwas_table=gwas_table.annotate(AF=mt.rows()[gwas_table.locus, gwas_table.alleles].variant_QC_Hail.AF[1])
        gwas_table=gwas_table.annotate(AC=mt.rows()[gwas_table.locus, gwas_table.alleles].pheno_call_stats.AC[1])
        gwas_table=gwas_table.annotate(AN=mt.rows()[gwas_table.locus, gwas_table.alleles].pheno_call_stats.AN)
        gwas_table=gwas_table.annotate(AF=mt.rows()[gwas_table.locus, gwas_table.alleles].pheno_call_stats.AF[1])
        gwas_table=gwas_table.annotate(n_HomRef=mt.rows()[gwas_table.locus, gwas_table.alleles].pheno_call_stats.homozygote_count[0])
        gwas_table=gwas_table.annotate(n_Het=mt.rows()[gwas_table.locus, gwas_table.alleles].pheno_n_het)
        gwas_table=gwas_table.annotate(n_HomAlt=mt.rows()[gwas_table.locus, gwas_table.alleles].pheno_call_stats.homozygote_count[1])
        gwas_table=gwas_table.key_by('locus')
        gwas_table=gwas_table.drop('alleles')
        gwas_table=gwas_table.select('rsid','REF','ALT','n', 'AF', 'beta', 'standard_error', 'p_value','nmr_phenotypes', 'n_HomRef','n_Het', 'n_HomAlt' )
        print(" Writing gwas table checkpoint")
        
        gwas = gwas_table.checkpoint(f"{tmp_dir}/gwas/{project}-{dataset}-gwas-nmr-{nmr2_new[0]}.table", overwrite=True)
        
        print("Exporting tsv table")
        gwas.export(f"{tmp_dir}/gwas/{project}-{dataset}-gwas-nmr-{nmr2_new[0]}.tsv.bgz", header=True)
        for j in range(len(nmr_new)):
            print(f"Plotting manhattan {j}:{nmr2_new[j]}")
            p = hl.plot.manhattan(gwas.p_value, title=f"{nmr2_new[j]} GWAS")
            output_file(f"{temp_dir}/gwas/{project}-{dataset}-{nmr2_new[j]}-manhattan.html", mode='inline')
            save(p)
            print(f"Plotting QQ plot for {j} - {nmr2_new[j]}")    
            q = hl.plot.qq(gwas.p_value, collect_all=False, n_divisions=100, title=f"{nmr2_new[j]} QQ plot")
            output_file(f"{temp_dir}/gwas/{project}-{dataset}-{nmr2_new[j]}-QQplot.html", mode='inline')
            save(q)
        nmr_new=[]
        nmr2_new=[]
            
        

        
            # print(f"Plotting manhattan plot for {index} - {pheno_name}")
            # p = hl.plot.manhattan(gwas.p_value, title=f"{pheno_name} GWAS")
            # output_file(f"{temp_dir}/gwas/{project}-{dataset}-{index}-{i}-manhattanplot.html", mode='inline')
            # save(p)
            # print(f"Plotting QQ plot for {index} - {pheno_name}")    
            # q = hl.plot.qq(gwas.p_value, collect_all=False, n_divisions=100, title=f"{pheno_name} QQ plot")
            # output_file(f"{temp_dir}/gwas/{project}-{dataset}-{index}-{i}-QQplot.html", mode='inline')
            # save(q)
