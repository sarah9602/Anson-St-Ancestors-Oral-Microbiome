#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: sarah9602

Will get a single genome file in kraken format given a single taxid.

USAGE
    get_fasta_in_kraken_format.py /PATH/TO/INPUT.fa TAXID
"""

import os,sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

## Get directory information and input file
if len(sys.argv[1].split('/')) > 1:
    DIR='/'.join(sys.argv[1].split('/')[:-1])
    FILE=sys.argv[1].split('/')[-1]
else:
    DIR=os.getcwd()
    FILE=sys.argv[1]

if not sys.argv[2]:
    print("Please provide a taxid")
    sys.exit(1)
else:
    TAX=sys.argv[2]

os.chdir(DIR)

def get_fasta_in_kraken_format(inputfile,taxid):
    # output=open(outfile,'a+')
    ID="_".join(inputfile.split("_")[:2])
    print("Converting %s to kraken format" % ID)
    outfile='.'.join(inputfile.split('.')[:-1])+'.tax.fa'
    o=open(outfile,'a+')
    for seq_record in SeqIO.parse(inputfile,"fasta"):
        ID=seq_record.id
        NEWID="%s|kraken:taxid|%s" % (ID,taxid)
        newseq=SeqRecord(seq=seq_record.seq,id=NEWID,description=seq_record.description)
        SeqIO.write(newseq,o,'fasta')
    o.close()
    # Optional - remove input file
    # os.remove(inputfile)
    return

get_fasta_in_kraken_format(FILE,TAX)
