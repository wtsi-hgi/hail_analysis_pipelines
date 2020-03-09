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

def annotate_allele_counts_chrX(mt):
    #annotate samples with gender information from sampleQC table, change "ID" and "SEX"column names appropriately:
    #mt=matrixtable.annotate_cols(sex=sample_info_table.key_by("ID")[matrixtable.s].SEX)
    #calculate call stats(AC,AN, AF filtered by sex and if it's in PAR region)
    #mt=mt.annotate_rows(call_stats_male= hl.agg.filter(mt.sex == 1, hl.agg.call_stats(mt.GT, mt.alleles)))
    Number_males = mt.aggregate_cols(hl.agg.count_where(mt.sex ==1 ))
    Number_females = mt.aggregate_cols(hl.agg.count_where(mt.sex ==2 ))
    mt=mt.annotate_rows(call_stats_male_NONPAR= hl.agg.filter(((mt.sex == 1) & (mt.locus.in_x_nonpar())), hl.agg.call_stats(mt.GT, mt.alleles)))
    mt=mt.annotate_rows(call_stats_male_PAR= hl.agg.filter(((mt.sex == 1) & (mt.locus.in_x_par())), hl.agg.call_stats(mt.GT, mt.alleles)))
    mt=mt.annotate_rows(call_stats_female= hl.agg.filter(mt.sex == 2, hl.agg.call_stats(mt.GT, mt.alleles)))
    #Halve AC for males and halve AN for males and create 2 new fields maleAC, maleAN
    mt=mt.annotate_rows(maleAC_NONPAR = mt.call_stats_male_NONPAR.AC.map(lambda x:hl.int32(x/2)))
    mt=mt.annotate_rows(maleAN_NONPAR = hl.int32(mt.call_stats_male_NONPAR.AN / 2))
    #Add the new AC and AN for males with the AC,AN for females together to get the final count per variant:
    #final_AC and final_AN have the final values
    ## Calculating total allele number, which will be required for call rate calculation
    mt=mt.annotate_rows(Total_allele=hl.cond((mt.locus.in_x_nonpar()), ((Number_females*2)+ Number_males), ((Number_females + Number_males)*2)))
    
    mt=mt.annotate_rows(final_AC= mt.maleAC_NONPAR + mt.call_stats_male_PAR.AC + mt.call_stats_female.AC)
    mt=mt.annotate_rows(final_AN= mt.maleAN_NONPAR + mt.call_stats_male_PAR.AN + mt.call_stats_female.AN) 
    mt=mt.annotate_rows(frequency=hl.float64(mt.final_AC[1]/mt.final_AN))
    mt=mt.annotate_rows(final_call_rate=hl.float64((mt.maleAN_NONPAR + mt.call_stats_male_PAR.AN + mt.call_stats_female.AN)/(mt.Total_allele)))

    fields_to_drop=['call_stats_male_NONPAR','call_stats_male_PAR','call_stats_female','maleAC_NONPAR','maleAN_NONPAR', 'Total_allele']
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

    #Read chrX matrixtable before final variant QC 
    #(2894876, 11863)
    mtX=hl.read_matrix_table(f"{temp_dir}/checkpoints/chrX_sampleQC_step1_filtered_GenotypeQC_checkpoint.mt")
    mtX=annotate_allele_counts_chrX(mtX)
    mtX_qc = hl.variant_qc(mtX, name='variant_QC_Hail')
    mtX_qc_filtered = mtX_qc.filter_rows(
        (mtX_qc.final_call_rate >= 0.95) &
        (mtX_qc.variant_QC_Hail.gq_stats.mean >= 20) &
        (mtX_qc.final_AC[1] >= 1)
    )
    # mtX_qc_filtered.count() : (2705791, 11863)
    mtX_qc_filtered=mtX_qc_filtered.checkpoint(f"{tmp_dir}/intervalwgs/chrX_final_AC_corrected_filtered_variant_QC.mt", overwrite=True)
    #2020-03-09 15:53:06 Hail: INFO: wrote matrix table with 2705791 rows and 11863 columns in 900 partitions 
    #to hdfs://spark-master:9820//intervalwgs/chrX_final_AC_corrected_filtered_variant_QC.mt