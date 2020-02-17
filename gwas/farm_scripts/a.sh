#!/bin/bash 
conda init bash
conda activate pa10
python3 /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/scripts/hail_analysis_pipelines/gwas/explode_gwas_table.py --table /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/big-am.gz --number am

#!/bin/bash 
conda init bash
conda activate pa10
python3 /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/scripts/hail_analysis_pipelines/gwas/explode_gwas_table.py --table /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/big-an.gz --number an

#!/bin/bash 
conda init bash
conda activate pa10
python3 /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/scripts/hail_analysis_pipelines/gwas/explode_gwas_table.py --table /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/big-ao.gz --number ao


#!/bin/bash 
conda init bash
conda activate pa10
python3 /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/scripts/hail_analysis_pipelines/gwas/explode_gwas_table.py --table /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/big-ap.gz --number ap

#!/bin/bash 
conda init bash
conda activate pa10
python3 /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/scripts/hail_analysis_pipelines/gwas/explode_gwas_table.py --table /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/big-aq.gz --number aq


#!/bin/bash 
conda init bash
conda activate pa10
python3 /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/scripts/hail_analysis_pipelines/gwas/explode_gwas_table.py --table /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/big-ar.gz --number ar

#!/bin/bash 
conda init bash
conda activate pa10
python3 /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/scripts/hail_analysis_pipelines/gwas/explode_gwas_table.py --table /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/big-as.gz --number as

#!/bin/bash 
conda init bash
conda activate pa10
python3 /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/scripts/hail_analysis_pipelines/gwas/explode_gwas_table.py --table /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/big-at.gz --number at

#!/bin/bash 
conda init bash
conda activate pa10
python3 /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/scripts/hail_analysis_pipelines/gwas/explode_gwas_table.py --table /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/big-aw.gz --number aw

#!/bin/bash 
conda init bash
conda activate pa10
python3 /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/scripts/hail_analysis_pipelines/gwas/explode_gwas_table.py --table /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/big-ax.gz --number ax


#!/bin/bash 
conda init bash
conda activate pa10
python3 /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/scripts/hail_analysis_pipelines/gwas/explode_gwas_table.py --table /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/big-ay.gz --number ay


#!/bin/bash 
conda init bash
conda activate pa10
python3 /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/scripts/hail_analysis_pipelines/gwas/explode_gwas_table.py --table /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/big-az.gz --number az


bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutaa.o -e joberroraa.e ./explode.sh aa
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutab.o -e joberrorab.e ./explode.sh ab
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutac.o -e joberrorac.e ./explode.sh ac 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutad.o -e joberrorad.e ./explode.sh ad 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutae.o -e joberrorae.e ./explode.sh ae 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutaf.o -e joberroraf.e ./explode.sh af 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutag.o -e joberrorag.e ./explode.sh ag 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutah.o -e joberrorah.e ./explode.sh ah 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutai.o -e joberrorai.e ./explode.sh ai 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutaj.o -e joberroraj.e ./explode.sh aj 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutak.o -e joberrorak.e ./explode.sh ak 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutal.o -e joberroral.e ./explode.sh al 

bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutam.o -e joberroram.e ./explode.sh am 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutan.o -e joberroran.e ./explode.sh an 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutao.o -e joberrorao.e ./explode.sh ao 

bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutap.o -e joberrorap.e ./explode.sh ap 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutaq.o -e joberroraq.e ./explode.sh aq 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutar.o -e joberrorar.e ./explode.sh ar 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutas.o -e joberroras.e ./explode.sh as 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutat.o -e joberrorat.e ./explode.sh at 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutau.o -e joberrorau.e ./explode.sh au 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutav.o -e joberrorav.e ./explode.sh av 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutaw.o -e joberroraw.e ./explode.sh aw 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutax.o -e joberrorax.e ./explode.sh ae 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutay.o -e joberroray.e ./explode.sh ay 
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutaz.o -e joberroraz.e ./explode.sh az




