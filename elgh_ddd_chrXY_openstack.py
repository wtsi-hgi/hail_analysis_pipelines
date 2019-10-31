
'''
author: Pavlos Antoniou
date: 22/07/19
'''

import os
import hail as hl
import pyspark
import json
import sys




project_root = os.path.dirname(os.path.dirname(__file__))
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

    #####################################################################
    ###################### chromosome X  ##############################
    #####################################################################
    #DEFINE INPUT FILE PATHS
    chrX_haploid_vcf_path= storage["elghddd"]["s3"]["chrXhap"]
    chrX_diploid_vcf_path= storage["elghddd"]["s3"]["chrXdip"]
    chrY_haploid_vcf_path= storage["elghddd"]["s3"]["chrYhap"]
    chrY_diploid_vcf_path= storage["elghddd"]["s3"]["chrYdip"]
    print("Read in all VCFs:")
    mtX_hap = hl.import_vcf(chrX_haploid_vcf_path, force_bgz=True, reference_genome='GRCh38')
    mtX_dip = hl.import_vcf(chrX_diploid_vcf_path, force_bgz=True, reference_genome='GRCh38')
    mtY_hap = hl.import_vcf(chrY_haploid_vcf_path, force_bgz=True, reference_genome='GRCh38')
    mtY_dip = hl.import_vcf(chrY_diploid_vcf_path, force_bgz=True, reference_genome='GRCh38')
    ############################
    CHROMOSOME = "chrX"
    print(f"Reading {CHROMOSOME} diploid VCF calls")
    mtX_dip_split = hl.split_multi_hts(mtX_dip, keep_star = False)
    mtX_hap_split = hl.split_multi_hts(mtX_hap, keep_star = False)
    ####### sex imputation
    print("Sex imputation:")
    mtx_unphased = mtX_dip_split.select_entries(GT=hl.unphased_diploid_gt_index_call(mtX_dip_split.GT.n_alt_alleles()))
    imputed_sex = hl.impute_sex(mtx_unphased.GT)

    #Haploid VCF: filter to male samples only
    print("Haploid chrX VCF: filter to male samples only")
    mtX_hap_males = mtX_hap_split.filter_cols(imputed_sex[mtX_hap_split.s].is_female, keep=False)
    #Diploid vcf: one matrixtable only Female, one matrixtable only male
    print("Diploid chrX vcf: create new matrixtable one for males one for female samples")
    mtX_dip_males = mtX_dip_split.filter_cols(imputed_sex[mtX_dip_split.s].is_female, keep = False)
    mtX_dip_females = mtX_dip_split.filter_cols(imputed_sex[mtX_dip_split.s].is_female, keep = True)

    #Get par region
    rg = hl.get_reference('GRCh38')
    par = rg.par

    #For haploid male chrX variants remove all falling within par
    mtX_hap_males_NONPAR = hl.filter_intervals(mtX_hap_males, par, keep=False)

    #For diploid male samples chrX variants KEEP PAR regions only
    mtX_dip_males_PAR = hl.filter_intervals(mtX_dip_males, par, keep = True)

    #The hap and dip VCF files differ in entry fielse. Dropping these fields to allow union_rows and union_cols to proceed
    fields_to_drop = ['PGT', 'PID', 'PS']

    mtX_dip_males_PAR_dropf = mtX_dip_males_PAR.drop(*fields_to_drop)
    mtX_dip_females_dropf = mtX_dip_females.drop(*fields_to_drop)
    mtX_union_males = mtX_hap_males_NONPAR.union_rows(mtX_dip_males_PAR_dropf)


    mt_final = mtX_union_males.union_cols(mtX_dip_females_dropf)
    print("Output matrixtable:")
    mt_final = mt_final.checkpoint(
        f"{tmp_dir}/elgh-ddd/elgh-ddd-chrX_surgery.mt", overwrite=True)
    print("Export VCF")
    hl.export_vcf(mt_final, f"{tmp_dir}/elgh-ddd/elgh-ddd-chrX-surgery.vcf.bgz")

