
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

#read tsv file in pandasl

#'["blbal"]'.strip("[]")

tsv1="/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/10lines_long.txt.gz"
tsvout="/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/test.txt.bgz"
def explode(df, lst_cols, fill_value=''):
    # make sure `lst_cols` is a list
    if lst_cols and not isinstance(lst_cols, list):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)

    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()

    if (lens > 0).all():
        # ALL lists in cells aren't empty
        return pd.DataFrame({
            col:np.repeat(df[col].values, df[lst_cols[0]].str.len())
            for col in idx_cols
        }).assign(**{col:np.concatenate(df[col].values) for col in lst_cols}) \
          .loc[:, df.columns]
    else:
        # at least one list in cells is empty
        return pd.DataFrame({
            col:np.repeat(df[col].values, df[lst_cols[0]].str.len())
            for col in idx_cols
        }).assign(**{col:np.concatenate(df[col].values) for col in lst_cols}) \
          .append(df.loc[lens==0, idx_cols]).fillna(fill_value) \
          .loc[:, df.columns]

if __name__ == "__main__":

    chunk_list=[]
    df_chunk = pd.read_csv(tsv1, delimiter="\t",compression='gzip', chunksize=5)

    for chunk in df_chunk:
    #df= pd.read_csv(tsv1, delimiter="\t", compression='gzip')
        chunk=chunk.drop(columns=['alleles'])
       # chunk[['beta','standard_error','p_value','nmr_phenotypes']]=chunk[['beta','standard_error','p_value','nmr_phenotypes']].apply(lambda x: x.str.strip("[]"), axis=1)
       #explode(chunk, lst_cols=['beta','standard_error','p_value','nmr_phenotypes'])
        chunk=chunk.set_index(['locus','rsid','n','REF','ALT','AF']).apply(lambda x: x.apply(pd.Series).stack()).reset_index(level=1, drop=True)
        #chunk=chunk.explode(['beta','standard_error','p_value','nmr_phenotypes'])
       # chunk[['nmr_phenotypes']]=chunk[['nmr_phenotypes']].apply(lambda x: x.str.strip("\"\""), axis=1)
        #chunk_list.append(chunk)
        #chunk.columns=['locus',	'alleles',	'rsid',	'n', 'beta','standard_error','p_value',	'nmr_phenotypes','REF',	'ALT','AF']
        
        chunk.to_csv(tsvout, sep="\t", header=(not os.path.exists(tsvout)),compression='gzip', index=False ,mode='a' )

    #df_concat=pd.concat(chunk_list)
    #df.to_csv(tsvout, sep="\t",compression='gzip', header=True, index=False )







