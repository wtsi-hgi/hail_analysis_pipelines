
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



def fix_genotpe(mt):

    pl_expr = (
        hl.zip_with_index(mt.PL).flatmap(lambda x: hl.range(x[0] + 1).map(lambda i: hl.cond(i == x[0], x[1], 999))))

    gt_allele = mt.GT.n_alt_alleles()
    gt_and_pl = (hl.switch(mt.GT.ploidy).when(1, (hl.call(gt_allele, gt_allele), pl_expr)).default((mt.GT, mt.PL)))

    mt = mt.annotate_entries(GT=gt_and_pl[0], PL=gt_and_pl[1])

    return mt

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
    #chrY_diploid_vcf_path= storage["elghddd"]["s3"]["chrYdip"]
    print("Read in all VCFs:")
    mtX_hap = hl.import_vcf(chrX_haploid_vcf_path, force_bgz=True, reference_genome='GRCh38')
    mtX_dip = hl.import_vcf(chrX_diploid_vcf_path, force_bgz=True, reference_genome='GRCh38')
    mtY_hap = hl.import_vcf(chrY_haploid_vcf_path, force_bgz=True, reference_genome='GRCh38')
    #mtY_dip = hl.import_vcf(chrY_diploid_vcf_path, force_bgz=True, reference_genome='GRCh38')
    ############################
    #Haploid VCF fixing of GT
    print("Fix GT haploid")
    mtX_hap=fix_genotpe(mtX_hap)
    mtY_hap=fix_genotpe(mtY_hap)
    mtX_hap = mtX_hap.checkpoint(f"{tmp_dir}/elgh-ddd/chrX_haploid_GT_fixed.mt", overwrite=True)
    mtY_hap = mtY_hap.checkpoint(f"{tmp_dir}/elgh-ddd/chrY_haploid_GT_fixed.mt", overwrite=True)






    ############################### chrX
    CHROMOSOME = "chrX"
    print(f"Split multi_hts {CHROMOSOME} diploid VCF calls")
    mtX_dip_split = hl.split_multi_hts(mtX_dip, keep_star = False)
    mtX_hap_split = hl.split_multi_hts(mtX_hap, keep_star = False)
    mtX_dip_split = mtX_dip_split.checkpoint(f"{tmp_dir}/elgh-ddd/chrX_diploid_split.mt", overwrite=True)
    mtX_hap_split = mtX_hap_split.checkpoint(f"{tmp_dir}/elgh-ddd/chrX_haploid_split.mt", overwrite=True)



    ####### sex imputation
    print("Sex imputation:")
    mtx_unphased = mtX_dip_split.select_entries(GT=hl.unphased_diploid_gt_index_call(mtX_dip_split.GT.n_alt_alleles()))
    imputed_sex = hl.impute_sex(mtx_unphased.GT)

    ################## chrY
    print("chrY  surgery")
    mtY_hap_split = hl.split_multi_hts(mtY_hap, keep_star=False)
    call1 = hl.eval(hl.call(0, 0))
    mtY_hap_sex_annot = mtY_hap_split.annotate_cols(
        sex=hl.cond(imputed_sex[mtY_hap_split.s].is_female, "female", "male"))
    mtY_hap_sex_annot = mtY_hap_sex_annot.annotate_entries(
        GT=hl.cond(mtY_hap_sex_annot.sex == "female", call1, mtY_hap_sex_annot.GT))
    print("Checkpoint chrY")
    mtY_hap_sex_annot= mtY_hap_sex_annot.checkpoint(f"{tmp_dir}/elgh-ddd/chrY_haploid_final.mt", overwrite=True)
    print("Export chrY VCF")
    hl.export_vcf(mtY_hap_sex_annot, f"{tmp_dir}/elgh-ddd/elgh-ddd-chrY-surgery_final.vcf.bgz")

    #Haploid VCF: filter to male samples only
    print("Haploid chrX VCF: filter to male samples only")
    mtX_hap_males = mtX_hap_split.filter_cols(imputed_sex[mtX_hap_split.s].is_female, keep=False)
    mtX_hap_males = mtX_hap_males.checkpoint(f"{tmp_dir}/elgh-ddd/chrX_haploid_males.mt", overwrite=True)

    #Diploid vcf: one matrixtable only Female, one matrixtable only male
    print("Diploid chrX vcf: create new matrixtable one for males one for female samples")
    mtX_dip_males = mtX_dip_split.filter_cols(imputed_sex[mtX_dip_split.s].is_female, keep = False)
    mtX_dip_females = mtX_dip_split.filter_cols(imputed_sex[mtX_dip_split.s].is_female, keep = True)
    mtX_dip_males = mtX_dip_males.checkpoint(f"{tmp_dir}/elgh-ddd/chrX_diploid_males.mt", overwrite=True)
    mtX_dip_females = mtX_dip_females.checkpoint(f"{tmp_dir}/elgh-ddd/chrX_diploid_females.mt", overwrite=True)


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
    hl.export_vcf(mt_final, f"{tmp_dir}/elgh-ddd/elgh-ddd-chrX-surgery_final.vcf.bgz")

