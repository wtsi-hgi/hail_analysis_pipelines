
'''
author: Pavlos Antoniou
date: 22/07/19
'''

import os
import hail as hl
import pyspark
import json
import sys




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

    #####################################################################
    ###################### chromosome X  ##############################
    #####################################################################
    #DEFINE INPUT FILE PATHS

    gwas1=hl.read_table(f"{tmp_dir}/intervalwgs/WGS-autosomes-gwasfbc-checkpoint-somalogic_0.05_partitioned")


    gwas_table=gwas1.flatten()

    gwas_table.export(f"{tmp_dir}/intervalwgs/gwas-wgsautosomes-export-somalogic_p_value_0.05.tsv.bgz", header=True)