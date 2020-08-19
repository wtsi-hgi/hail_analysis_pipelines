#!/bin/bash 

# bsub -q long -G hgi -R'select[mem>10000] rusage[mem=10000]' -M10000  -o gwas_bed_addpheno.o -e gwas_bed_addpheno.e  ./gwas_to_bed_addpheno.sh 

folder=/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/final_results/gwas_tables/nmr
outdir=${folder}/bed_files
mkdir -p $outdir

for f in ${folder}/*.tsv.bgz
do 
fbname=$(basename "$f" .tsv.bgz)
zless $f  |tr ':' '\t' > ${outdir}/interim.bed
awk 'BEGIN{FS=OFS="\t"} {$2 = $2 OFS $2} 1' ${outdir}/interim.bed| tail -n+2| gzip -c > ${outdir}/${fbname}.bed.gz
rm ${outdir}/interim.bed
done