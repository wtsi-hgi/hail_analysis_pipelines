import pandas as pd 
import os
import json
import pathlib
from pathlib import Path
import re
#open path 
plots_path="/lustre/scratch119/humgen/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/final_tables"
allpaths=[p for p in pathlib.Path(plots_path).iterdir() if p.is_file()]

#print(allpaths)

data=[]


for path in allpaths:
    print(path)
    info={}
    filename=path.stem
    m=re.search(r'INT-WGS-gwas-nmr-(.*).(\w+).tsv.gz', filename)
    if m:
        name=m.group(1)
        print(name)
        
    
    #print(filename)

