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

def annotate_allele_counts_chrY(mt):
    #annotate samples with gender information from sampleQC table, change "ID" and "SEX"column names appropriately:
    #mt=matrixtable.annotate_cols(sex=sample_info_table.key_by("ID")[matrixtable.s].SEX)
    #calculate call stats(AC,AN, AF filtered by sex and if it's in PAR region)
    #mt=mt.annotate_rows(call_stats_male= hl.agg.filter(mt.sex == 1, hl.agg.call_stats(mt.GT, mt.alleles)))
    Number_males = mt.aggregate_cols(hl.agg.count_where(mt.sex ==1 ))
    Number_females = mt.aggregate_cols(hl.agg.count_where(mt.sex ==2 ))
    mt=mt.annotate_rows(call_stats_male= hl.agg.filter(mt.sex == 1 , hl.agg.call_stats(mt.GT, mt.alleles)))
    mt=mt.annotate_rows(call_stats_female= hl.agg.filter(mt.sex == 2, hl.agg.call_stats(mt.GT, mt.alleles)))
    #Halve AC for males and halve AN for males and create 2 new fields maleAC, maleAN
    mt=mt.annotate_rows(maleAC = mt.call_stats_male.AC.map(lambda x:hl.int32(x/2)))
    mt=mt.annotate_rows(maleAN = hl.int32(mt.call_stats_male.AN / 2))
    #Add the new AC and AN for males with the AC,AN for females together to get the final count per variant:
    #final_AC and final_AN have the final values
    mt=mt.annotate_rows(final_AC= mt.maleAC + mt.call_stats_female.AC)
    mt=mt.annotate_rows(final_AN= mt.maleAN + mt.call_stats_female.AN) 
    mt=mt.annotate_rows(frequency=hl.float64(mt.final_AC[1]/mt.final_AN))
    mt=mt.annotate_rows(final_call_rate=hl.float64((mt.maleAN + mt.call_stats_female.AN)/Number_males))
    fields_to_drop=['call_stats_male','call_stats_female','maleAC','maleAN']
    mt=mt.drop(*fields_to_drop)
    return mt;


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

    #Read chrY matrixtable before final variant QC 
    #(44390, 11863)
    mtY=hl.read_matrix_table(f"{temp_dir}/intervalwgs/sex_chromosomes/chrY_sampleQC_step1_filtered_GenotypeQC_checkpoint.mt")
    
    print("Variant qc for male samples:")
    mtY_males=mtY.filter_entries(mtY.sex==1,keep=True)
    mtY_females=mtY.filter_entries(mtY.sex==1,keep=False)
    mtY_males=mtY_males.filter_cols(mtY_males.sex==1,keep=True)
    mtY_males = hl.variant_qc(mtY_males, name='variant_QC_Hail')
    mtY_males_filtered = mtY_males.filter_rows(
        (mtY_males.variant_QC_Hail.call_rate >= 0.95) &
        (mtY_males.variant_QC_Hail.gq_stats.mean >= 20) &
        (mtY_males.variant_QC_Hail.AC[1] >= 1)
    )
    mt_sqc_vqc_filtered=mtY_males_filtered.union_cols(mtY_females)
    #####################################################################
    ###################### FINAL QC AFTER FILTERING  ####################
    #####################################################################

    fields_to_drop = ['variant_QC_Hail', 'sample_QC_Hail']

    mt1 = mt_sqc_vqc_filtered.drop(*fields_to_drop)

    mt2 = hl.sample_qc(mt1, name='sample_QC_Hail')
    mt3 = hl.variant_qc(mt2, name='variant_QC_Hail')



    mtY_corrected=annotate_allele_counts_chrY(mt3)
    
    # mtX_qc_filtered.count() : (10515, 11863)
    mtY_corrected=mtY_corrected.checkpoint(f"{tmp_dir}/intervalwgs/chrY_final_AC_corrected_filtered_variant_QC.mt", overwrite=True)
    #2020-03-09 16:42:31 Hail: INFO: wrote matrix table with 10515 rows and 11863 columns in 250 partitions to hdfs://spark-master:9820//intervalwgs/chrY_final_AC_corrected_filtered_variant_QC.mt
    