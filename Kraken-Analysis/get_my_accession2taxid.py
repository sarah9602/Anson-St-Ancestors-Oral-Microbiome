#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: sarah9602 (2024)

Create a mapping file with NCBI accession, accession version, taxids, and GI number. Extracted from larger nucl_wgs.accession2taxid.gz acquired from NCBI. 

Usage
    get_my_accession2taxid.py list_of_accession_numbers existing_mapping_file

List of accession numbers can be obtained from:
    awk '{print $1}' plastid.name.map > accnums.txt;awk '{print $1}' mitochondrion.name.map >> accnums.txt
Existing mapping file can be nucl_wgs.accession2taxid.gz from NCBI ftp.
To reduce file size, grep "_" from this file and output condensed.accession2taxid.
"""
import os,sys
import pandas as pd

## Get directory information and input file
if len(sys.argv[1].split('/')) > 1:
    DIR='/'.join(sys.argv[1].split('/')[:-1])
    FILE=sys.argv[1].split('/')[-1]
else:
    DIR=os.getcwd()
    FILE=sys.argv[1]
    
os.chdir(DIR)
acclis=[x.rstrip() for x in open(FILE,'r').readlines()]
## Huge file - lots of space
acc2taxid=pd.read_csv(sys.argv[2],sep='\t',index_col=False)
acc2taxidhere=acc2taxid[acc2taxid["accession.version"].isin(acclis)]
with open('my.accession2taxid.map','w') as o:
    acc2taxidhere.to_csv(o,sep='\t',index=False)
    