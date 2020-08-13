
'''
author: Pavlos Antoniou
#chrX non PAR regions to convert to haploid for male samples before exporting to VCF. 
# For chrY convert all variants to haploid before and export VCF.
date: 06/08/20
'''

import os
import hail as hl
import pyspark
import json
import sys
from pathlib import Path


project_root = Path(__file__).parent.parent.parent
print(project_root)

s3credentials = os.path.join(project_root, "config_files/s3_credentials.json")
print(s3credentials)

storage = os.path.join(project_root, "config_files/storage.json")

thresholds = os.path.join(project_root, "config_files/thresholds.json")

with open(f"{s3credentials}", 'r') as f:
    credentials = json.load(f)

with open(f"{storage}", 'r') as f:
    storage = json.load(f)

with open(f"{thresholds}", 'r') as f:
    thresholds = json.load(f)


def fix_genotype(mt):

    pl_expr = (
        hl.zip_with_index(mt.PL).flatmap(lambda x: hl.range(x[0] + 1).map(lambda i: hl.cond(i == x[0], x[1], 999))))

    gt_allele = mt.GT.n_alt_alleles()
    gt_and_pl = (hl.switch(mt.GT.ploidy).when(
        1, (hl.call(gt_allele, gt_allele), pl_expr)).default((mt.GT, mt.PL)))

    mt = mt.annotate_entries(GT=gt_and_pl[0], PL=gt_and_pl[1])

    return mt


if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory
    tmp_dir = "hdfs://spark-master:9820/"
    temp_dir = os.path.join(os.environ["HAIL_HOME"], "tmp")
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])

    #####################################################################
    ###################### chromosome X  ##############################
    #####################################################################
    # chrX non PAR regions to convert to haploid for male samples before exporting to VCF.
    # For chrY convert all variants to haploid before and export VCF.
    # Read chromosome mt
    # select male samples
    # select non par regions and convert GT to haploid

    mtX_dip = hl.read_matrix_table(
        f"{temp_dir}/intervalwgs/sex_chromosomes/chrX_final_AC_corrected_filtered_variant_QC_ordered_FINAL.mt")
    mt = mtX_dip
    mt = mt.annotate_entries(GT=hl.case()
                             .when((mt.GT.is_diploid())
                                   & (mt.sex == 1)
                                   & (mt.locus.in_x_nonpar()), hl.call(mt.GT[0], phased=False))
                             .default(hl.call(mt.GT[0], mt.GT[1], phased=False)))
    mt.write(f"{temp_dir}/intervalwgs/sex_chromosomes/chrX_final_AC_corrected_filtered_variant_QC_ordered_FINAL_haploid_males_non_par.mt")
    hl.export_vcf(
        mt, f"{tmp_dir}/VCFs/chrX_final_AC_corrected_filtered_variant_QC_ordered_FINAL_haploid_males_non_par.vcf.bgz")

    #####################################################################
    ###################### chromosome Y     ##############################
    #####################################################################
    mtY_dip = hl.read_matrix_table(
        f"{temp_dir}/intervalwgs/sex_chromosomes/chrY_final_AC_corrected_filtered_variant_QC_ordered_FINAL.mt")

    mt = mtY_dip
    mt = mt.annotate_entries(GT=hl.call(mt.GT[0], phased=False))
    mt.write(f"{temp_dir}/intervalwgs/sex_chromosomes/chrY_final_AC_corrected_filtered_variant_QC_ordered_FINAL_haploid.mt")
    hl.export_vcf(
        mt, f"{tmp_dir}/VCFs/chrY_final_AC_corrected_filtered_variant_QC_ordered_FINAL_haploid.vcf.bgz")
