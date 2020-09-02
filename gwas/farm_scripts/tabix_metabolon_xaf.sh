#!/bin/bash 

# bsub -q long -G hgi -R'select[mem>10000] rusage[mem=10000]' -M10000  -o gwas_bed_addpheno.o -e gwas_bed_addpheno.e  ./gwas_to_bed_addpheno.sh 
# 
CHUNK=/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/final_results/gwas_tables/metabolon/bed_files/xaf

folder=/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/final_results/gwas_tables/metabolon/bed_files
outdir=${folder}
mkdir -p $outdir
while read -r FILE; do
	COUNT=$(( $COUNT + 1 ))
	fbname=$(basename "${FILE}" .bed.gz)
    gzip -c -d $FILE | bgzip -c > ${outdir}/${fbname}.bed.bgz
    tabix -pbed ${outdir}/${fbname}.bed.bgz
	echo ${COUNT}
done < "${CHUNK}"

