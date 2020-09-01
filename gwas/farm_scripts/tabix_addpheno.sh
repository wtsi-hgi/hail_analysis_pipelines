#!/bin/bash 

# bsub -q long -G hgi -R'select[mem>10000] rusage[mem=10000]' -M10000  -o gwas_bed_addpheno.o -e gwas_bed_addpheno.e  ./gwas_to_bed_addpheno.sh 
# mercury 
folder=/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/final_results/gwas_tables/addpheno/bed_files
outdir=${folder}/bed_files
mkdir -p $outdir

for f in ${folder}/*.bed.gz
do 
tabix -p $f 
done