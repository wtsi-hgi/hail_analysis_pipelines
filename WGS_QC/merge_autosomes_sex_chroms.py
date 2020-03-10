'''
Coorect AC and AN for chrX and chrY and run QC
author: Pavlos Antoniou
date: 09/03/2020
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

    dbsnp = hl.import_vcf("s3a://intervalwgs-qc/dbsnp_153.hg38.vcf.gz", force_bgz=True, reference_genome='GRCh38',
                          skip_invalid_loci=True)
    dbsnp_rows = dbsnp.rows()
    mt = hl.read_matrix_table(f"{temp_dir}/intervalwgs/WGS_final_february_2020_updated_rsID.mt")
    #Remove chrX and chrY
    mt=mt.filter_rows(mt.locus.contig == "chrX",keep=False)
    mt=mt.filter_rows(mt.locus.contig == "chrY",keep=False)
    #Read new mtX and mtY
    mtX=hl.read_matrix_table(f"{temp_dir}/intervalwgs/sex_chromosomes/chrX_final_AC_corrected_filtered_variant_QC.mt")
    mtY=hl.read_matrix_table(f"{temp_dir}/intervalwgs/sex_chromosomes/chrY_final_AC_corrected_filtered_variant_QC.mt")
    mt=mt.drop(mt.call_stats)
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
    info1 = mtY.info.drop(*info_fields_to_drop)
    mtY = mtY.annotate_rows(info=info1)
    print("Reorder rows to match mt")
    #Reorder rows to match mt
    mt=mt.add_row_index()
    mtX=mtX.add_row_index()
    mtX=mtX.select_rows(*mt.row_value)
    mtY=mtY.add_row_index()
    mtY=mtY.select_rows(*mt.row_value)

    #Reorder cols to match mt
    print("Reorder cols to match mt")
    samples=mt.s.collect()
    name_to_index = {col['s']: col['col_idx'] for col in
    mtX.add_col_index().cols().key_by().select('col_idx', 's').collect()}
    newmtX = mtX.choose_cols([name_to_index[s] for s in samples])
    name_to_index = {col['s']: col['col_idx'] for col in
    mtY.add_col_index().cols().key_by().select('col_idx', 's').collect()}
    newmtY = mtY.choose_cols([name_to_index[s] for s in samples])

    #Join 
    print("Join mts")
    mt=mt.union_rows(newmtX)
    mt=mt.union_rows(newmtY)
    print("Annotating dbsnp")
    mt=mt.annotate_rows(rsid=dbsnp_rows[mt.locus, mt.alleles].rsid)

    #Write to disk
    print("writing chrY ")
    newmtY=newmtY.checkpoint(f"{tmp_dir}/intervalwgs/sex_chromosomes/chrY_final_AC_corrected_filtered_variant_QC_ordered_FINAL.mt")
    print("writing chrX ")
    newmtX=newmtX.checkpoint(f"{tmp_dir}/intervalwgs/sex_chromosomes/chrX_final_AC_corrected_filtered_variant_QC_ordered_FINAL.mt")
    print("writing WGS ")
    mt=mt.checkpoint(f"{tmp_dir}/intervalwgs/WGS_final_march_2020_dbsnp_v53.mt")