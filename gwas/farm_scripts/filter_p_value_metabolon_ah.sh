#!/bin/bash 
# bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o fiilter_pvalue.o -e filter_pvalue.e  


CHUNK=/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/final_results/gwas_tables/metabolon/xah_original
folder=/lustre/scratch119/humgen/projects/interval_wgs/analysis/hail_analysis/gwas/final_results/gwas_tables/metabolon
outfolder=${folder}/pvalue_filter
mkdir -p ${outfolder}
while read -r FILE; do
	COUNT=$(( $COUNT + 1 ))
	fbname=$(basename "${FILE}" .tsv.bgz)
    zcat ${FILE} |  awk 'BEGIN {FS="\t";OFS="\t"} {if (NR ==1 || $8 < 0.00000005)  {print}}' | gzip > ${outfolder}/${fbname}.p5e-8.tsv.bgz; 
	echo ${COUNT}
done < "${CHUNK}"