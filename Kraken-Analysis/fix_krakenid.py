#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Replace sequence headers in kraken-formatted genome files reflecting new TAXIDs created using "taxonkit create-taxdump".
Usage:
    fix_krakenid.py INFILE /PATH/TO/NEW/TAXDUMP/taxid.map

@author: sarah9602 (2024)
"""

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

FASTA=sys.argv[1]
TAXMAP=sys.argv[2]

taxidmap=dict(x.rstrip().split(None,1) for x in open(TAXMAP,'r'))

outfile='.'.join(FASTA.split(".fa")[:-1]) +'.formatted.fa'

def fix_kraken_id(inputfile,outfile="genomes.fa"):
    # output=open(outfile,'a+')
    DIR='/'.join(inputfile.split('/')[-1])
    o=open(outfile,'a+')
    x=open('%s/problem_taxids.txt' % DIR,'a+')
    for seq_record in SeqIO.parse(inputfile,"fasta"):
        oldtax=seq_record.id.split('|')[-1]
        try:
            newtax=taxidmap[oldtax]
        except:
            x.write('%s\n' % oldtax)
            continue
        NEWID=seq_record.id.replace(oldtax,newtax)
        desc=' '.join(seq_record.description.split(' ')[1:])
        newseq=SeqRecord(seq=seq_record.seq,id=NEWID,description=desc)
        SeqIO.write(newseq,o,'fasta')
    o.close()
    x.close()
    # os.remove(inputfile)
    return

fix_kraken_id(FASTA,outfile)


