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

d={}
i=1
for path in allpaths:
    print(path)
    df = pd.read_csv(path, delimiter="\t", nrows=2, compression='gzip')

    rows,columns=df.shape
    if columns==11:
        d[path]=columns


for k, v in d.items():
    print(k)
