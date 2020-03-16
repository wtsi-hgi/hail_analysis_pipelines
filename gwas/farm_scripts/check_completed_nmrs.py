import pandas as pd 
import os
import json
import pathlib
from pathlib import Path
import re
#open path 
plots_path="/lustre/scratch119/humgen/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/gwas_tables"
allpaths=[p for p in pathlib.Path(plots_path).iterdir() if p.is_file()]

#print(allpaths)

data=[]
nmr_phenotypes='/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/gwas_tables/nmr_pheno_names.txt'
with open(nmr_phenotypes) as f:
    lines = [line.rstrip() for line in f]

for path in allpaths:
    #print(path)
    info={}

    filename=path.stem
    print(path.stem)
    m=re.search(r'INT-WGS-gwas-nmr-(.*).tsv.bgz', filename)
    if m:
        name=m.group(1)
        #print(name)
        if name in data:
           data.append(name)

print(data)
for ph in lines:
    if ph not in data:
        print(ph)

    #print(filename)

