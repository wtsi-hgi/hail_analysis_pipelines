
'''
QC for 1 chromosome
author: Pavlos Antoniou
date: 22/07/19
'''

import os
import hail as hl
import sys
import json



project_root = os.path.dirname(os.path.dirname(__file__))
print(project_root)



snp_vqsr_threshold= -0.6647
indel_vqsr_threshold=-1.2537

def fix_genotpe(mt):

    pl_expr = (
        hl.zip_with_index(mt.PL).flatmap(lambda x: hl.range(x[0] + 1).map(lambda i: hl.cond(i == x[0], x[1], 999))))

    gt_allele = mt.GT_original.n_alt_alleles()
    gt_and_pl = (hl.switch(mt.GT_original.ploidy).when(1, (hl.call(gt_allele, gt_allele), pl_expr)).default((mt.GT_original, mt.PL)))

    mt = mt.annotate_entries(GT=gt_and_pl[0], PL=gt_and_pl[1])

    return mt

BUCKET = "gs://interval-wgs"
#Define chromosome here
tmp_dir="/Users/pa10/Programming/google-code/google/tmp"

if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    #Define the hail persistent storage directory
    hl.init(default_reference="GRCh38", tmp_dir=tmp_dir)



    #DEFINE INPUT FILE PATHS
    chrX=hl.read_matrix_table('gs://interval-wgs/chrX.mt')
    chrY=hl.read_matrix_table('gs://interval-wgs/chrY.mt')

    #mtY_dip = hl.import_vcf(chrY_diploid_vcf_path, force_bgz=True, reference_genome='GRCh38')
    ############################
    #Haploid VCF fixing of GT
    print("renaming original GT")
    chrX_renamed=chrX.rename({'GT':'GT_original'})
    chrY_renamed = chrY.rename({'GT': 'GT_original'})
    print("Fix GT haploid")
    #chrX=fix_genotpe(chrX_renamed)
    chrY=fix_genotpe(chrY_renamed)
    print("Writin out matrixtables")
    CHROMOSOME="chrX"
    chrX.write(f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}_GT_fixed.mt", overwrite=True)
    CHROMOSOME = "chrY"
    chrY.write(f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}_GT_fixed.mt",
               overwrite=True)
