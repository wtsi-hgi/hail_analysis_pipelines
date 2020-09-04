#!/bin/bash 

# bsub -q long -G hgi -R'select[mem>10000] rusage[mem=10000]' -M10000  -o gwas_bed_addpheno.o -e gwas_bed_addpheno.e  ./gwas_to_bed_addpheno.sh 
# 
CHUNK=/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/final_results/gwas_tables/nmr/bed_files/failed_tbis.fofn

folder=/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/final_results/gwas_tables/nmr/bed_files
outdir=${folder}
mkdir -p $outdir
while read -r FILE; do
	COUNT=$(( $COUNT + 1 ))
	fbname=$(basename "${FILE}" _corrected.bed.bgz)
    sortBed -i  ${FILE}  > ${outdir}/${fbname}_sorted.bed
    bgzip -c ${outdir}/${fbname}_sorted.bed > ${outdir}/${fbname}_sorted.bed.bgz
    tabix -pbed ${outdir}/${fbname}_sorted.bed.bgz
	echo ${COUNT}
done < "${CHUNK}"

