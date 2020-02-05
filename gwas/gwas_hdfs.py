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
     
    phenotypes = hl.import_table("s3a://intervalwgs-qc/phenotypes_pc_fbc_nmr_olink_somalogic_16-12-1019.csv", impute=True, delimiter=',')
    dbsnp = hl.import_vcf("s3a://intervalwgs-qc/GCF_000001405.38.fixed.vcf.gz", force_bgz=True, reference_genome='GRCh38',skip_invalid_loci=True)
    dbsnp_rows = dbsnp.rows()

    CHROMOSOME="WGS"
    mt = hl.read_matrix_table(f"{temp_dir}/intervalwgs/WGS_final_january_2020_updated.mt")
    
    print("Number of initial variants:")
    print(mt.count())
    
    print("Annotating matrixtable with phenotypes:")
    phenotypes=phenotypes.key_by('ID')
    mt = mt.annotate_cols(phenotype=phenotypes[mt.s])

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


    print("Linear regression")
    for index,value in enumerate(nmr):
        print(nmr2[index])
        pheno_name=nmr2[index]
        gwas = hl.linear_regression_rows(
            y=value,
            x=mt.GT.n_alt_alleles(), covariates=[1.0]+pcas[0:10], pass_through=[mt.rsid])
        
        
        #Filter p-value
        #gwas1=gwas_annotated.filter(gwas_annotated.p_value < 5e-8 , keep=True)
        # gwas1 = gwas1.checkpoint(f"{tmp_dir}/gwas/gwas{pheno_name}-pvalue5e-8.table", overwrite=True)
        print(" Writing gwas table checkpoint")
        gwas = gwas.checkpoint(f"{tmp_dir}/gwas/{project}-{dataset}-gwas-{pheno_name}.table", overwrite=True)
        # print(gwas_annotated.count())
        #print(gwas1.count())
        #gwas1=gwas.filter(gwas.p_value[0].any(lambda x: x < 5e-8 ), keep=True)
        #gwas1.export(f"{tmp_dir}/gwas/gwas-{pheno_name}_pvalue-5e-8.tsv.bgz", header=True)
        #ANNOTATE WITH DBSNP
        gwas_annotated = gwas.annotate(dbsnp=dbsnp_rows[gwas.locus, gwas.alleles].rsid)
        print("Exporting tsv table")
        gwas.export(f"{tmp_dir}/gwas/{project}-{dataset}-gwas-{pheno_name}.tsv.bgz", header=True)
       # gwas_table = hl.import_table(f"{tmp_dir}/gwas/gwas-{pheno_name}_test_loop_pvalue-5e-8.tsv.bgz", key=['locus', 'alleles'],
          #                       types={'locus': 'locus<GRCh38>', 'alleles': 'array<str>'})
          #                       types={'locus': 'locus<GRCh38>', 'alleles': 'array<str>'})
        
        print(f"Plotting manhattan plot for {index} - {pheno_name}")
        p = hl.plot.manhattan(gwas.p_value, title=f"{pheno_name} GWAS")
        output_file(f"{temp_dir}/gwas/{project}-{dataset}-{index}-{pheno_name}-manhattanplot.html", mode='inline')
        save(p)
        print(f"Plotting QQ plot for {index} - {pheno_name}")    
        q = hl.plot.qq(gwas.p_value, collect_all=False, n_divisions=100, title=f"{pheno_name} QQ plot")
        output_file(f"{temp_dir}/gwas/{project}-{dataset}-{index}-{pheno_name}-QQplot.html", mode='inline')
        save(q)
