
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
import gzip
import numpy as np
import ast 
from ast import literal_eval
import argparse


#read tsv file in pandasl

#'["blbal"]'.strip("[]")
parser = argparse.ArgumentParser()

parser.add_argument("--table", type=str, help="path to the gwas table to be exploded")
parser.add_argument("--number", type=str, help="id for outfile")

args=parser.parse_args()
outpath="/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/final_tables"
id_out=args.number
if __name__ == "__main__":

    chunk_list=[]
    #locus   alleles rsid    n       beta    standard_error  p_value nmr_phenotypes  REF     ALT     AF
    #locus   rsid    REF     ALT     n       AF      beta    standard_error  p_value nmr_phenotypes
    colnames=['locus','rsid','REF','ALT','n','AF','beta','standard_error','p_value','nmr_phenotypes']
    df_chunk = pd.read_csv(args.table, delimiter="\t", names=colnames, compression='gzip', chunksize=500000, header=None, low_memory=False)
    
    for chunk in df_chunk:
        #print(chunk['rsid'])
        #print(chunk)
        #print("before")
    #df= pd.read_csv(tsv1, delimiter="\t", compression='gzip')
        chunk['beta'] = chunk['beta'].apply(literal_eval)
        chunk['standard_error'] = chunk['standard_error'].apply(literal_eval)
        chunk['p_value'] = chunk['p_value'].apply(literal_eval)
        chunk['nmr_phenotypes'] = chunk['nmr_phenotypes'].apply(literal_eval)
        #chunk=chunk.drop(columns=['alleles'])
        #print("after:")
        #print(chunk)
        
       # chunk[['beta','standard_error','p_value','nmr_phenotypes']]=chunk[['beta','standard_error','p_value','nmr_phenotypes']].apply(lambda x: x.str.strip("[]"), axis=1)
        #explode(chunk, lst_cols=['beta','standard_error','p_value','nmr_phenotypes'])
        chunk=chunk.set_index(['locus','rsid' ,'REF','ALT','n','AF']).apply(lambda x: x.apply(pd.Series).stack()).reset_index()
        #locus, rsid, REF, ALT, n, AF, beta, standard_error, p_value, nmr_phenotypes
        #chunk=chunk.drop(columns=['level_5'])
        #chunk=chunk.reset_index()
        #print(chunk)
        #chunk=chunk.explode(['beta','standard_error','p_value','nmr_phenotypes'])
       # chunk[['nmr_phenotypes']]=chunk[['nmr_phenotypes']].apply(lambda x: x.str.strip("\"\""), axis=1)
        #chunk_list.append(chunk)
        #chunk.columns=['locus',	'alleles',	'rsid',	'n', 'beta','standard_error','p_value',	'nmr_phenotypes','REF',	'ALT','AF']
        
        for i, x in chunk.groupby('nmr_phenotypes'):
            p = os.path.join(outpath, "INT-WGS-gwas-nmr-{}.{}.tsv.gz".format(i.lower(), id_out))
            x.to_csv(p, sep="\t", header=(not os.path.exists(p)),compression='gzip',index=False, mode='a')

        #chunk.to_csv(tsvout, sep="\t", header=(not os.path.exists(tsvout)),compression='gzip', index=True ,mode='a' )

    #df_concat=pd.concat(chunk_list)
    #df.to_csv(tsvout, sep="\t",compression='gzip', header=True, index=False )







