#!/bin/bash

conda init bash
conda activate pa10


file=$1
name=$2
zcat ${file} | awk '!seen[$0]++' | gzip > /lustre/scratch119/humgen/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/gwas_tables/remove_duplicates/${name}