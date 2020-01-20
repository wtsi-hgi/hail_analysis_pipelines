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

phenotypes = hl.import_table("s3a://intervalwgs-qc/phenotypes_pc_fbc_nmr_olink_somalogic_16-12-1019.csv", impute=True, delimiter=',')
dbsnp = hl.import_vcf("s3a://intervalwgs-qc/GCF_000001405.38.fixed.vcf.gz", force_bgz=True, reference_genome='GRCh38',
                          skip_invalid_loci=True)
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
     
    CHROMOSOME="WGS-autosomes"
    mt = hl.read_matrix_table(f"{temp_dir}/matrixtables/WGS_autosomes-full-sampleqc-variantqc-FILTERED.mt")
    
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


    print("Linear regression")
    for index,value in enumerate(nmr[0:10]):
        print(nmr2[index])
        pheno_name=nmr2[index]
        gwas = hl.linear_regression_rows(
            y=value,
            x=mt.GT.n_alt_alleles(), covariates=[1.0]+pcas[0:10], pass_through=[mt.rsid])

        gwas_annotated = gwas_table.annotate(dbsnp=dbsnp_rows[gwas_table.locus, gwas_table.alleles].rsid)
        #gwas1=gwas.filter(gwas.p_value[0].any(lambda x: x < 5e-8 ), keep=True)
        gwas1=gwas.filter(gwas.p_value < 5e-8 , keep=True)
       
        gwas1 = gwas1.checkpoint(f"{tmp_dir}/gwas/gwas{pheno_name}-test_pvalue5e-8.table", overwrite=True)
        print(gwas1.count())
        #gwas1=gwas.filter(gwas.p_value[0].any(lambda x: x < 5e-8 ), keep=True)
        gwas1.export(f"{tmp_dir}/gwas/gwas-{pheno_name}_test_loop_pvalue-5e-8.tsv.bgz", header=True)
        
        gwas_table = hl.import_table(f"{tmp_dir}/gwas/gwas-{pheno_name}_test_loop_pvalue-5e-8.tsv.bgz", key=['locus', 'alleles'],
                                 types={'locus': 'locus<GRCh38>', 'alleles': 'array<str>'})

        dbsnp_rows = dbsnp.rows()
        gwas_annotated = gwas_table.annotate(dbsnp=dbsnp_rows[gwas_table.locus, gwas_table.alleles].rsid)

    print("export table")
    gwas_annotated.export(f"{BUCKET}/gwas/gwas_WGS_with_rsID_not_parallel.tsv.bgz", delimiter="\t",header=True)
