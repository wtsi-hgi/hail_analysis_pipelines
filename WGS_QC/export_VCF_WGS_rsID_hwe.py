
'''
Export VCF from WGS interval
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

CHROMOSOMES= [
              # "chr1",
              # "chr2",
              # "chr3",
              # "chr4",
              # "chr5",
               "chr6",
               "chr7",
               "chr8",
               "chr9",
               "chr10",
               "chr11",
               "chr12",
               "chr13",
               "chr14",
               "chr15",
               "chr16",
               "chr17"

               #"chr18", #done
               #"chr19",
               #"chr20",
               #"chr21",
               #"chr22",
               #"chrX",
               #"chrY"
]

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
    
    mt_wgs=hl.read_matrix_table(f"{temp_dir}/intervalwgs/WGS_after_surgery-full-sampleqc-variantqc-FILTERED.mt")
   
    dbsnp_vcf= "s3a://intervalwgs-qc/dbsnp_153.hg38.vcf.gz"
    dbsnp = hl.import_vcf(dbsnp_vcf, force_bgz=True, reference_genome='GRCh38',
                          skip_invalid_loci=True)
    dbsnp_rows = dbsnp.rows()
    for CHROMOSOME in CHROMOSOMES:
        print(CHROMOSOME)
       
        mt=mt_wgs.filter_rows(mt_wgs.locus.contig==CHROMOSOME)
        mt = hl.variant_qc(mt)
        mt = mt.annotate_rows(info = mt.info.annotate(AC=mt.variant_qc.AC[1]))
        mt = mt.annotate_rows(info = mt.info.annotate(AF=mt.variant_qc.AF[1]))
        mt = mt.annotate_rows(info = mt.info.annotate(AN=mt.variant_qc.AN))
        mt = mt.annotate_rows(info= mt.info.annotate(pvalueHWE=mt.variant_qc.p_value_hwe))
        #mt = mt.annotate_rows(info=mt.info.annotate(dbsnp_v3=dbsnp_rows[mt.locus, mt.alleles].rsid))
        mt=mt.annotate_rows(rsid=dbsnp_rows[mt.locus, mt.alleles].rsid)
        #Export everything 
        #mt1 = mt.select_entries()
        hl.export_vcf(mt, f"{tmp_dir}/VCFs/{CHROMOSOME}/{CHROMOSOME}.intervalwgs_v1_all_info.vcf.bgz", parallel="separate_header")

        #Export only Genotype info in exonic regions
        mt2=mt.select_entries(mt.GT)
        #interval_table = hl.import_locus_intervals(chr_exon_regions, reference_genome='GRCh38')
        #filtered_mt = mt2.filter_rows(hl.is_defined(interval_table[mt2.locus]))
        hl.export_vcf(mt2, f"{tmp_dir}/VCFs/{CHROMOSOME}/{CHROMOSOME}.intervalwgs_v2_GT_only.vcf.bgz", parallel="separate_header")
    #Export everything but genotype info
    #hl.export_vcf(mt, f"{BUCKET}/VCFs/{CHROMOSOME}/{CHROMOSOME}.noGT_after_sampleQC.vcf.bgz")

    #fields_to_drop = ['variant_QC_Hail', 'sample_QC_Hail']

    #mt1 = mt.drop(*fields_to_drop)
    #mt2 = hl.sample_qc(mt1, name='sample_QC_Hail')
    #mt3 = hl.variant_qc(mt2, name='variant_QC_Hail')


    #mt3 = mt3.checkpoint(
    #    f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-full-sampleqc-variantqc-filtered-FINAL.mt", overwrite=True)
    #mt3_cols = mt3.cols()
    #mt3_cols.flatten().export(
    #    f"{BUCKET}/output-tables/{CHROMOSOME}/{CHROMOSOME}-sampleQC_filtered_FINAL.tsv.bgz", header=True)

    #mt3_rows = mt3.rows()
    #mt3_rows.select(mt3_rows.variant_QC_Hail).flatten().export(
    #    f"{BUCKET}/output-tables/{CHROMOSOME}/{CHROMOSOME}-variantQC_filtered_FINAL.tsv.bgz", header=True)

