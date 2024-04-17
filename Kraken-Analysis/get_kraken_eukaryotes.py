#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 15:35:05 2024

Get top eukaryotic signatures from merged krkaen2mpa file
Usage:
    get_kraken_eukaryotes.py INFILE

@author: sarahjane
"""
import sys,os
import pandas as pd

## Get directory information and input file
if len(sys.argv[1].split('/')) > 1:
    DIR='/'.join(sys.argv[1].split('/')[:-1])
    FILE=sys.argv[1].split('/')[-1]
else:
    DIR=os.getcwd()
    FILE=sys.argv[1]

os.chdir(DIR)
df=pd.read_csv(FILE,sep='\t',index_col=False)
df['average']=df.mean(numeric_only=True,axis=1)
dfsort=df.sort_values(by=['average'],ascending=False)
dftop=dfsort.head(n=20).drop(columns=['average'])
with open('%s/top_20_euk.txt' % DIR,'w') as o:
    for i in list(dftop['#Classification'].str.split('s__').str[-1]):
        o.write('%s\n' % i)
outfile="%s.top20.txt" % '.'.join(FILE.split('.')[:-1])
with open(outfile,'w') as x:
    dftop.to_csv(x,sep='\t',index=False)