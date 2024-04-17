#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 15:32:57 2024

Inspect sourcetracker feature tables
Output will be put in input directory

@author: sarahjane
"""
import os,sys,glob,argparse
import pandas as pd


argp = argparse.ArgumentParser( prog = "sourcetracker_feature_contributions.py",
    description = """Output tables with feature information from sourcetracker2.""")
argp.add_argument( "-s",    metavar = "/PATH/TO/STUDY/sourcetracker2.RARE.taxonlevel", nargs = 1,
    help = "Path to sourcetracker2 output directory." )
argp.add_argument( '-st',    metavar = "/PATH/TO/SOURCETRACKER/METAPHLAN", nargs = 1,
    help = "Path to SourceTracker Metaphlan output directory." )
argp.add_argument( '-t',    metavar = "taxonlevel", nargs = 1,
    help = "Taxonomic level sourcetracker was run" )

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" )

argp.usage = argp.format_usage()[7:]+"""\n\n\tFirst run sourcetracker2 on Metaphlan4 output files with --per_sink_feature_assignments flag."""
# DIR="/Users/sarahjane/Documents/PhD/Isabella/sourcetracker2.30000.p"
# os.chdir(DIR)
args = argp.parse_args( )
otumap=dict(x.rstrip().split(None,1) for x in open('%s/otu-map-%s.txt' % (args.st[0],args.t[0].lower()),'r'))

def filter_table(table):
    out=table.loc[table[table.iloc[:,1:] > 5].any(axis=1)]
    return out

def format_table(infile):
    df=pd.read_csv(infile,sep='\t',index_col=0)
    dft=df.transpose()
    dft['taxonomy']=dft.soindex.map(otumap)
    cols=['taxonomy']+list(dft.columns)[:-1]
    dfout=dft[cols]
    out=filter_table(dfout)
    return out
    
for i in glob.glob('%s/*feature_table.txt' % args.s[0]):
    samp=i.split('/')[-1].split('.')[0]
    df=format_table(i)
    with open('%s/%s.top_contributors.mapped.txt' % (args.s[0],samp), 'w') as x:
        df.to_csv(x,sep='\t',index=False)
