
'''
format gwas output on farm5
author: Pavlos Antoniou
date: 10/02/2020
conda activate pa10
'''

import os
import json
import sys
from pathlib import Path
import pandas as pd 


#read tsv file in pandasl

#'["blbal"]'.strip("[]")

tsv1="/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/INT-WGS-gwas-nmr-1.tsv.bgz"
tsvout="/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/INT-WGS-gwas-nmr-acace-final.tsv.bgz"
if __name__ == "__main__":

    #reader = pd.read_csv(file1, delimiter="\t",chunksize=128)
    #for chunk in reader:
    df= pd.read_csv(tsv1, delimiter="\t", compression='gzip')

    df=df.apply(lambda x: x.strip("[]"), axis=1)

    df.to_csv(tsvout, sep="\t",compression='gzip', header=True, index=False )

    





