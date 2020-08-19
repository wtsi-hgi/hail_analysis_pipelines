#!/bin/bash 

# bsub -q long -G hgi -R'select[mem>10000] rusage[mem=10000]' -M10000  -o gwas_bed_addpheno.o -e gwas_bed_addpheno.e  ./gwas_to_bed_addpheno.sh 

CHUNK=/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/final_results/gwas_tables/metabolon/xah

folder=/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/final_results/gwas_tables/metabolon
outdir=${folder}/bed_files
mkdir -p $outdir
bed=interimxah.bed
while read -r FILE; do
	COUNT=$(( $COUNT + 1 ))
	fbname=$(basename "${FILE}" .tsv.bgz)
    zless "${FILE}"  |tr ':' '\t' > "${outdir}"/${bed}
	#s3cmd put "${FILE}" "${BUCKET}/${CHROMOSOME}-vcf/${gs_name}"
    awk 'BEGIN{FS=OFS="\t"} {$2 = $2 OFS $2} 1' ${outdir}/${bed} | tail -n+2| gzip -c > ${outdir}/${fbname}.bed.gz
    rm ${outdir}/${bed}
	echo ${COUNT}
done < "${CHUNK}"

