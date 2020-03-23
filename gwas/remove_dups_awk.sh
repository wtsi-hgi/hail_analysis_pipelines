#!/bin/bash

conda init bash
conda activate pa10



zcat /lustre/scratch119/humgen/projects/interval_wgs/analysis/ha
il_analysis/gwas/nmr_results/gwas_tables/INT-WGS-gwas-nmr-nmr_lldlpl_.tsv.bgz | awk '!seen[$0]++' | gzip > /lustre/scratch119/humgen/projects/interval_wgs/analysis/ha
il_analysis/gwas/nmr_results/gwas_tables/remove_duplicates/INT-WGS-gwas-nmr-nmr_lldlpl_.tsv.bgz 