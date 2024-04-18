#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make sure gtdb.taxid.map is in /PATH/TO/GTDB/DIRECTORIES
USAGE:
    python get_gtdb_in_kraken_format.py /PATH/TO/GTDB/DIRECTORIES

@author: sarah9602 (2024)
"""

# you'll see this alias in documentation, examples, etc.
import sys,subprocess,os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
 
# Get input directory or quit
if len(sys.argv[1].split('/')) > 1:
    DIR=sys.argv[1]
elif len(sys.argv[1].split('/')) == 1:
    DIR='%s/%s' % (os.getcwd(),sys.argv[1])
elif sys.argv[1] == ".":
    DIR=os.getcwd()
else:
    print('Provide a path to gtdb directory')
    print('Directory must contain gtdb.name.map and gtdb.taxid.map.')
    sys.exit()
os.chdir(DIR)

# Define function to alter header ids
def get_fasta_in_kraken_format(inputfile,TAX,outfile="genomes.fa"):
    # output=open(outfile,'a+')
    o=open(outfile,'a+')
    for seq_record in SeqIO.parse(inputfile,"fasta"):
        ID=seq_record.id
        NEWID="%s|kraken:taxid|%s" % (ID,TAX)
        newseq=SeqRecord(seq=seq_record.seq,id=NEWID,description=seq_record.description)
        SeqIO.write(newseq,o,'fasta')
    o.close()
    # os.remove(inputfile)
    return

## Make sure taxid map has GTDB taxids
gtdbtaxmap=dict(x.rstrip().split(None,1) for x in open('gtdb.taxid.map','r'))

# Walk through input directories and collect all fna files
for root,dirs,files in os.walk(DIR):
    for f in files:
        if f.startswith(('GCA','GCF')):
            if f.endswith("gz"):
                fil=os.path.join(root,f)
                subprocess.call("gunzip "+fil,shell=True)
                newfil='.'.join(fil.split('.')[:-1])
                # os.remove(fil)
            elif f.endswith('.fna'):
                newfil=os.path.join(root,f)
            filname=''.join(newfil.split('/')[-1].split('.fna')[:-1])
            taxid=gtdbtaxmap[filname]
            get_fasta_in_kraken_format(newfil,taxid,'gtdb.genomes.fa')

