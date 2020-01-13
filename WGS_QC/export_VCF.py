
'''
QC for WGS
author: Pavlos Antoniou
date: 22/07/19
'''

import os
import hail as hl
import sys
import json


CHROMOSOMES = ["chr2",
               "chr3",
               "chr4",
               "chr5",
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
               "chr17",
               "chr18",
               "chr19",
               "chr20",
               "chr21",
               "chr22",
      
               ]



BUCKET = "gs://interval-wgs"
#Define chromosome here
tmp_dir="/Users/pa10/Programming/google-code/google/tmp"

if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    #Define the hail persistent storage directory
    hl.init(default_reference="GRCh38", tmp_dir=tmp_dir)








    #mt_chr1 = hl.read_matrix_table(f"{BUCKET}/matrixtables/chr1/chr1-full-sampleqc-variantqc-filtered-FINAL.mt")

#    for CHROMOSOME in CHROMOSOMES:
 #       mt = hl.read_matrix_table(f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-full-sampleqc-variantqc-filtered-FINAL.mt")
  #      mt_chr1 = mt_chr1.union_rows(mt)


    for CHROMOSOME in CHROMOSOMES:
        chr_exon_regions=f"{BUCKET}/exonic-regions/{CHROMOSOME}.txt"
    
        mt = hl.read_matrix_table(f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-full-sampleqc-variantqc-filtered-FINAL.mt")
        mt = hl.variant_qc(mt)
        mt = mt.annotate_rows(info = mt.info.annotate(AC=mt.variant_qc.AC[1]))
        mt = mt.annotate_rows(info = mt.info.annotate(AF=mt.variant_qc.AF[1]))
        mt = mt.annotate_rows(info = mt.info.annotate(AN=mt.variant_qc.AN))
        #Export everything but genotype info
        mt1 = mt.select_entries()
        #hl.export_vcf(mt1, f"{BUCKET}/VCFs/{CHROMOSOME}/{CHROMOSOME}.noGT_after_sampleQC.vcf.bgz")

        #Export only Genotype info in exonic regions
        mt2=mt.select_entries(mt.GT)
        interval_table = hl.import_locus_intervals(chr_exon_regions, reference_genome='GRCh38')
        filtered_mt = mt2.filter_rows(hl.is_defined(interval_table[mt2.locus]))
        hl.export_vcf(filtered_mt, f"{BUCKET}/VCFs/{CHROMOSOME}/{CHROMOSOME}.GT_only_after_sampleQC_exonic.vcf.bgz")
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

