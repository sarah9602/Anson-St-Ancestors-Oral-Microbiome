#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
How to read count mapping file:
ls STDIR/*mpa|while read in; do \
    samp=`echo $in|awk -F"/" '{print $NF}'|awk -F"." '{print $1}'`; \
    ct=`grep "reads processed" $in|sed 's/#//'|awk '{print $1}'|bc`; \
    echo -e "$samp\t$ct"; done > ST.counts.txt

@author: sarahjane
"""
import argparse
import os
import pandas as pd
import numpy as np

## Initialize output file to have ST (sourcetracker) merged with samples
## Relative abundnace converted to count data
def get_outfile(samp):
    FILE=samp.split('/')[-1]
    OUTFILE='ST-%s.counts.%s' % ('.'.join(FILE.split('.')[:-1]),FILE.split('.')[-1])
    return OUTFILE

## Main function
def merge( samp, st, sampmap ):
    """
    Outputs the table join of the given pre-split string collection.
    :param  samp:       merged sample relative metaphlan file
    :type   samp:       input file
    :param  st:         merged count sourcetracker file.
    :type   st:         input file
    :param  sampmap:    sample count map.
    :type   sampmap:    input file
    """
    if len(samp.split('/')) > 1:
        DIR='/'.join(samp.split('/')[:-1])
        FILE=samp.split('/')[-1]
    else:
        DIR=os.getcwd()
        FILE=samp
    
    # List sample file header
    sampcols=list(pd.read_csv('%s/%s' % (DIR,FILE),nrows=1,sep='\t'))
    # Load merged metaphlan table into dataframe
    sampdf=pd.read_csv('%s/%s' % (DIR,FILE),sep='\t',comment="#",index_col=0,usecols=[i for i in sampcols if i != "clade_taxid"])
    # Load sourcetracker table into dataframe
    stdf=pd.read_csv(st,sep='\t',index_col=0)
    # Load sample count map into dictionary
    sampct=dict(x.rstrip().split(None,1) for x in open(sampmap,'r'))
    
    # Initialize new sample dataframe with count data
    newsamp=pd.DataFrame()
    for a, b in sampdf.items():
        tot=int(sampct[a])
        b2=b.apply(lambda x: round(x/100*tot))
        if newsamp.empty:
            newsamp=b2.to_frame()
        else:
            newsamp=pd.merge(newsamp,b2.to_frame(),how='outer',left_index=True,right_index=True)
    
    # Merge sourcetracker count table with sample count table
    bigdf=pd.merge(newsamp,stdf,how='outer',left_index=True,right_index=True).fillna(0).apply(np.int64)
    return bigdf


argp = argparse.ArgumentParser( prog = "merge_metaphlan_sourcetracker_counts.py",
    description = """Joins metaphlan count data for samples with count data sourcetracker.""")
argp.add_argument( "-s",    metavar = "merged.SAMPLE.txt", nargs = 1,
    help = "Sample Merged metaphlan relative abundance." )
argp.add_argument( '-st',    metavar = "merged.ST.counts.txt", nargs = 1,
    help = "Merged SourceTracker counts file. Can be created using get_sourcetracker_metaphlan_counts.py script." )
argp.add_argument( '-m',    metavar = "SAMPLE.counts.txt", nargs = 1,
    help = "Map of sample name to read count. Can be acquired with metaphlan output." )

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" )

argp.usage = argp.format_usage()[7:]+"""\n\n\tPlease make sure to supply file paths to the files to combine.\n\n\tCan create mapping file as follows:\n\n\tls STDIR/*mpa|while read in; do \n\tsamp=`echo $in|awk -F"/" '{print $NF}'|awk -F"." '{print $1}'`; \n\tct=`grep "reads processed" $in|sed 's/#//'|awk '{print $1}'|bc`; \n\techo -e "$samp\t$ct"; done > ST.counts.txt"""


def main( ):
    args = argp.parse_args( )
    df=merge(args.s[0],args.st[0],args.m[0])
    outfile=get_outfile(args.s[0])
    with open(outfile,'w') as o:
        df.to_csv(o,sep='\t',index=True)

if __name__ == '__main__':
    main()