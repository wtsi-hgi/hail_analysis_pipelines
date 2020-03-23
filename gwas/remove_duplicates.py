
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
parser.add_argument("--name", type=str, help="id for outfile")

args=parser.parse_args()
outpath="/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/hail_analysis/gwas/nmr_results/tables/final_tables"
id_out=args.name
if __name__ == "__main__":

    chunk_list=[]
    #locus   alleles rsid    n       beta    standard_error  p_value nmr_phenotypes  REF     ALT     AF
    #locus   rsid    REF     ALT     n       AF      beta    standard_error  p_value nmr_phenotypes
    
    df = pd.read_csv(args.table, delimiter="\t", names=colnames, compression='gzip',  low_memory=False)
    df=df.drop_duplicates()
    
    p = os.path.join(outpath, "INT-WGS-gwas-nmr-{}.removed_dups.tsv.gz".format(id_out))
    df.to_csv(p, sep="\t", header=(not os.path.exists(p)),compression='gzip',index=False, mode='a')

        #chunk.to_csv(tsvout, sep="\t", header=(not os.path.exists(tsvout)),compression='gzip', index=True ,mode='a' )

    #df_concat=pd.concat(chunk_list)
    #df.to_csv(tsvout, sep="\t",compression='gzip', header=True, index=False )







