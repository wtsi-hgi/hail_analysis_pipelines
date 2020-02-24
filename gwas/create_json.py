import pandas as pd 
import os
import json
import pathlib
from pathlib import Path
import re
#open path 
plots_path="/Users/pa10/Programming/gwas_display_app/public/data/img"
allpaths=[p for p in pathlib.Path(plots_path).iterdir() if p.is_file()]

#print(allpaths)

data=[]

i=1
for path in allpaths:
    #print(path)
    info={}
    filename=path.stem
    m=re.search(r'INT-(W[G|E]S)-(.*)-', filename)
    if m:
        
        project=m.group(1)
        print(project)
        name=m.group(2)
        print(name)
        
    if "QQplot" in filename:
        info['id']=i
        info['project']=project
        info['category']='nmr'
        info['name']=name
        qqplot="data/img/"+filename+".html"
        #print(qqplot)
        info['qqplot']=qqplot
        info['manhattan']="data/img/"+f"INT-{project}-{name}-manhattan.html"
        i=i+1
    data.append(info)
    #print(filename)


data=list(filter(None, data))
print(data)
with open('/Users/pa10/Programming/gwas_display_app/public/data/data.json', 'w') as fout:
    json.dump(data , fout)