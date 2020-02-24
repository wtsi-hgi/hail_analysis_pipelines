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


instance=1
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutaa.o -e joberroraa.e ./explode.sh aa $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutab.o -e joberrorab.e ./explode.sh ab $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutac.o -e joberrorac.e ./explode.sh ac  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutad.o -e joberrorad.e ./explode.sh ad  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutae.o -e joberrorae.e ./explode.sh ae  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutaf.o -e joberroraf.e ./explode.sh af  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutag.o -e joberrorag.e ./explode.sh ag  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutah.o -e joberrorah.e ./explode.sh ah  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutai.o -e joberrorai.e ./explode.sh ai  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutaj.o -e joberroraj.e ./explode.sh aj  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutak.o -e joberrorak.e ./explode.sh ak  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutal.o -e joberroral.e ./explode.sh al  $instance

bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutam.o -e joberroram.e ./explode.sh am  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutan.o -e joberroran.e ./explode.sh an  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutao.o -e joberrorao.e ./explode.sh ao  $instance

bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutap.o -e joberrorap.e ./explode.sh ap  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutaq.o -e joberroraq.e ./explode.sh aq  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutar.o -e joberrorar.e ./explode.sh ar  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutas.o -e joberroras.e ./explode.sh as  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutat.o -e joberrorat.e ./explode.sh at  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutau.o -e joberrorau.e ./explode.sh au  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutav.o -e joberrorav.e ./explode.sh av  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutaw.o -e joberroraw.e ./explode.sh aw  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutax.o -e joberrorax.e ./explode.sh ae  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutay.o -e joberroray.e ./explode.sh ay  $instance
bsub -q long -G hgi -R'select[mem>30000] rusage[mem=30000]' -M30000  -o joboutaz.o -e joberroraz.e ./explode.sh az $instance

