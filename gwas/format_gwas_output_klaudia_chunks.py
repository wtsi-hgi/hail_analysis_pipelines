
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
    df_chunk = pd.read_csv(tsv1, delimiter="\t",compression='gzip', chunksize=1000000)
    for chunk in df_chunk:
    #df= pd.read_csv(tsv1, delimiter="\t", compression='gzip')

        chunk=chunk.apply(lambda x: x.str.strip("[]"), axis=1)
        #chunk_list.append(chunk)
        chunk.columns=['locus',	'alleles',	'rsid',	'n', 'beta','standard_error','p_value',	'nmr_phenotypes','REF',	'ALT','AF']
        chunk.to_csv(tsvout, sep="\t",compression='gzip', header=False, index=False ,mode='a' )

    #df_concat=pd.concat(chunk_list)
    #df.to_csv(tsvout, sep="\t",compression='gzip', header=True, index=False )

    





