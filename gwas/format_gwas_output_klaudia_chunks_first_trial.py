
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

    chunk_list=[]
    df_chunk = pd.read_csv(tsv1, delimiter="\t",compression='gzip', chunksize=100000)

    for chunk in df_chunk:
    #df= pd.read_csv(tsv1, delimiter="\t", compression='gzip')
        chunk=chunk.drop(columns=['alleles'])
        chunk[['beta','standard_error','p_value','nmr_phenotypes']]=chunk[['beta','standard_error','p_value','nmr_phenotypes']].apply(lambda x: x.str.strip("[]"), axis=1)
        
        chunk[['nmr_phenotypes']]=chunk[['nmr_phenotypes']].apply(lambda x: x.str.strip("\"\""), axis=1)

        
        chunk.to_csv(tsvout, sep="\t", header=(not os.path.exists(tsvout)),compression='gzip', index=False ,mode='a' )
    





