
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

project_root=Path(__file__).parent.parent.parent



#project_root = os.path.dirname(os.path.dirname(__file__))
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

def import_vcf(vcf_path,partitions):
    mt = hl.import_vcf(vcf_path,
                            force_bgz=True,
                            reference_genome='GRCh38',
                            array_elements_required=False,
                            min_partitions=partitions)

    if mt.n_partitions() > partitions:
        mt = mt.naive_coalesce(partitions)
    return mt


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
    
    
    sample_QC_nonHail = hl.import_table(storage["intervalwgs"]["s3"]["sampleQC_non_hail"], impute=True)

    partitions=250
    chrYhapVCF= "s3a://intervalwgs/haploid/chrY"
    mt_hap=import_vcf(chrYhapVCF,partitions)
    mtY_hap= mtY_hap.checkpoint(f"{tmp_dir}/intervalwgs/haploidY.mt", overwrite=True)

    mtY_hap=fix_genotpe(mtY_hap)
    mtY_hap = mtY_hap.checkpoint(f"{tmp_dir}/intervalwes/chrY_haploid_GT_fixed.mt", overwrite=True)

    #####################################################################
    ###################### chromosome Y  ##############################
    #####################################################################    
    print("chrY  surgery")
    mtY_hap_split = hl.split_multi_hts(mtY_hap, keep_star=False)

    mtY_hap_sex_annot = mtY_hap_split.annotate_cols(sex=sample_QC_nonHail.key_by("ID")[mtY_hap_split.s].SEX)
    
    mtY_hap_sex_annot = mtY_hap_sex_annot.annotate_entries(
        GT=hl.or_missing(mtY_hap_sex_annot.sex == 1, mtY_hap_sex_annot.GT))
    


    info_fields_to_drop = ['AS_BaseQRankSum',
                           'AS_FS',
                           'AS_InbreedingCoeff',
                           'AS_MQ',
                           'AS_MQRankSum',
                           'AS_QD',
                           'AS_ReadPosRankSum',
                           'AS_SOR',
                           'RAW_MQandDP',
                           'DB']
    #info1 = mtY_hap_sex_annot.info.drop(*info_fields_to_drop)
    #mtY_hap_sex_annot = mtY_hap_sex_annot.annotate_rows(info=info1)
    print("Checkpoint chrY")
    mtY_hap_sex_annot= mtY_hap_sex_annot.checkpoint(f"{tmp_dir}/intervalwgs/chrY_haploid_final.mt", overwrite=True)
    print("Export chrY VCF")
    hl.export_vcf(mtY_hap_sex_annot, f"{tmp_dir}/intervalwgs/chrY-surgery_final.vcf.bgz")