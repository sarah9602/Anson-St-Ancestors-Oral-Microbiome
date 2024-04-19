#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: sarah9602 (2024)

Get merged metaphlan relative abundance SourceTracker table into count data.
provide as input merged.ST.txt and mapping file of sample names to read count
Output will be input file name with 'counts' added before file suffix.
Output will also have taxid lineage column removed.

Requires pandas installed

How to read count mapping file:
ls STDIR/*mpa|while read in; do \
    samp=`echo $in|awk -F"/" '{print $NF}'|awk -F"." '{print $1}'`; \
    ct=`grep "reads processed" $in|sed 's/#//'|awk '{print $1}'|bc`; \
    echo -e "$samp\t$ct"; done > ST.counts.txt
    
USAGE
    get_sourcetracker_metaphlan_counts.py /PATH/TO/MERGED_ST.TXT /PATH/TO/MAPPING/FILE

"""
import pandas as pd
import sys,os

## Get directory information and input file
if len(sys.argv[1].split('/')) > 1:
    DIR='/'.join(sys.argv[1].split('/')[:-1])
    FILE=sys.argv[1].split('/')[-1]
else:
    DIR=os.getcwd()
    FILE=sys.argv[1]

## Initiate output file name
OUTFILE='%s.counts.%s' % ('.'.join(FILE.split('.')[:-1]),FILE.split('.')[-1])

try:
    stct=dict(x.rstrip().split(None,1) for x in open(sys.argv[2],'r'))
except:
    print("Your mapping file either doesn't exist or isn't in the correct format.")
    print("Please provide a file with sample names in col1 and read counts in col2.")
    sys.exit(1)

## Get header list
stcols=list(pd.read_csv('%s/%s' % (DIR,FILE),comment="#",nrows=1,sep='\t'))
## Read merged metaphlan table into dataframe
stdf=pd.read_csv('%s/%s' % (DIR,FILE),sep='\t',comment="#",index_col=0,usecols=[i for i in stcols if i != "clade_taxid"])

## Initiate new dataframe for count data
newst=pd.DataFrame()
for a, b in stdf.items():
    tot=int(stct[a])
    b2=b.apply(lambda x: round(x/100*tot))
    if newst.empty:
        newst=b2.to_frame()
    else:
        newst=pd.merge(newst,b2.to_frame(),how='outer',left_index=True,right_index=True)

## Write output to file
with open('%s/%s' % (DIR,OUTFILE),'w') as o:
    newst.to_csv(o,sep='\t',index=True)


    
