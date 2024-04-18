#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
USAGE:
    python get_organelle_in_kraken_format.py /PATH/TO/ORGANELLE/FNA/FILES
    organelle.taxid.map file must be present in this directory

@author: sarah9602 (2024)
"""
import sys,subprocess,os,glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
 
# DIR="/Users/sarahjane/Documents/PhD/GTDB"
if len(sys.argv[1].split('/')) > 1:
    DIR=sys.argv[1]
elif len(sys.argv[1].split('/')) == 1:
    DIR='%s/%s' % (os.getcwd(),sys.argv[1])
elif sys.argv[1] == ".":
    DIR=os.getcwd()
else:
    print('Provide a path to organelle directory')
    print('Directory must contain organelle.taxid.map')
    sys.exit()
os.chdir(DIR)
# FIL="plastid.3.genomic.gbff"

def get_fasta_in_kraken_format(inputfile,outfile="genomes.fa"):
    # output=open(outfile,'a+')
    taxmap=dict(x.rstrip().split(None,1) for x in open('organelle.taxid.map','r'))
    o=open(outfile,'a+')
    for seq_record in SeqIO.parse(inputfile,"fasta"):
        ID=seq_record.id
        taxid=taxmap[ID]
        NEWID="%s|kraken:taxid|%s" % (ID,taxid)
        newseq=SeqRecord(seq=seq_record.seq,id=NEWID,description=seq_record.description)
        SeqIO.write(newseq,o,'fasta')
    o.close()
    # os.remove(inputfile)
    return


for f in glob.glob('*.fna'):
    if f.startswith(('mitochondrion','plastid')):
        if f.endswith("gz"):
            fil='.'.join(f.split('.')[:-1])
            subprocess.call("gunzip "+f,shell=True)
        else:
            fil=f
        get_fasta_in_kraken_format(fil,'organelle_genomes.fa')
