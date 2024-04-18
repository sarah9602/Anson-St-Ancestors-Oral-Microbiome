#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: sarah9602 (2024)
info and error functions adopted frnm Aitor Blanco (aitor.blancomiguez@unitn.it)
Influenced by sejalmodha

Run
    python3 update_kraken_continue.py -h 
for usage information
"""

import pandas as pd
import glob
import subprocess
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import time
import sys

def info(s, init_new_line=True, exit=False, exit_value=0, time_stamp=True):
    if init_new_line:
        sys.stdout.write('\n')
    if time_stamp:
        timestamp = time.ctime(int(time.time()))
        sys.stdout.write('{} '.format(timestamp))
    sys.stdout.write('{}'.format(s))
    sys.stdout.write('\n')
    sys.stdout.flush()

    if exit:
        sys.exit(exit_value)

def error(s, init_new_line=False, exit=False, exit_value=1, time_stamp=True):
    if init_new_line:
        sys.stderr.write('\n')
    if time_stamp:
        timestamp = time.ctime(int(time.time()))
        sys.stderr.write('{} '.format(timestamp))
        
    sys.stderr.write('[e] {}\n'.format(s))
    sys.stderr.flush()

    if exit:
        sys.exit(exit_value)

def read_params():
    p = argparse.ArgumentParser(description='')

    p.add_argument('-p', '--path_to_krakendb', type=str, default='.',
                   help="Path to where the kraken db will be built.\nIf it doesn't exist, it will be created.\nDefault is the working directory.")

    p.add_argument('-c', '--taxonomic_category', type=str, default=None,
                   help="Taxonomic domain to download (Required).\nOptions include:\n\tbacteria, fungi, protozoa, archaea, viral")
    
    p.add_argument('-a', '--assembly_level', nargs='+', type=str, default=None,
                   help='List of desired assembly completion levels (Optional). Options include:\n\tComplete Genome (default), Chromosome, Scaffold, Contig')

    p.add_argument('-r', '--refseq_category', type=str, default=None,
                   help="Level of assembly completion (Optional). Options include:\n\trepresentative genome, reference genome\nTake care when choosing.")


    return p.parse_args()


def check_params(args):
    if not args.taxonomic_category:
        error('-c (or --taxonomic_category) must be specified', exit=True,
            init_new_line=True)
    elif not os.path.exists("%s/%s" % (args.path_to_krakendb,args.taxonomic_category)):
        os.makedirs("%s/%s" % (args.path_to_krakendb,args.taxonomic_category),exist_ok=True)

def check_done():
    donelis=[]
    for f in glob.glob('*_genomic.tax.fna'):
        donelis.append("_".join(f.split("/")[-1].split("_")[:2]))
    return donelis

def process_url_file(inputurl):
    file_suffix='_genomic.fna.gz'
    ID=inputurl.split("/")[-1]+file_suffix
    url="%s/%s" % (inputurl,ID)
    info("Downloading %s" % inputurl.split("/")[-1])
    subprocess.call("wget "+url,shell=True)
    subprocess.call("gunzip "+ID,shell=True)
    return      

def get_fasta_in_kraken_format(inputfile,taxid,outputfile='sequences.fa'):
    ID="_".join(inputfile.split("_")[:2])
    info("Converting %s to kraken format" % ID)
    # outfile=inputfile.replace("_genomic.fna","_genomic.tax.fna")
    o=open(outputfile,'a+')
    for seq_record in SeqIO.parse(inputfile,"fasta"):
        ID=seq_record.id
        NEWID="%s|kraken:taxid|%s" % (ID,taxid)
        newseq=SeqRecord(seq=seq_record.seq,id=NEWID,description=seq_record.description)
        SeqIO.write(newseq,o,'fasta')
    o.close()
    os.remove(inputfile)
    return

def download_genomes_continue(path_to_krakendb,taxonomic_category,*assembly_level,**kwargs):
    category = kwargs.get('category', None)
    os.chdir("%s/%s" % (path_to_krakendb,args.taxonomic_category))
    if os.path.exists('assembly_summary.txt'):
        info("Removing assembly summary.")
        os.remove('assembly_summary.txt')
    info("Downloading %s assembly summary." % taxonomic_category)
    subprocess.call("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/%s/assembly_summary.txt" % taxonomic_category,shell=True)
    info("If genomes have already been downloaded, assembly summary will be reduced.")
    alldone=check_done()
    df=pd.read_csv("assembly_summary.txt",sep='\t',skiprows=[0],index_col=False)
    dfneed=df[~df.iloc[:,0].isin(alldone)].reset_index(drop=True)
    dffinal=pd.DataFrame(columns=dfneed.columns)
    if None not in assembly_level:
        for a in assembly_level:
            print(a[0])
            info("limiting download to %s." % a[0])
            dfx=dfneed[dfneed.assembly_level==a[0]].reset_index(drop=True)
            dffinal=pd.concat([dffinal,dfx]).drop_duplicates().reset_index(drop=True)
    # else:
    #     dffinal=dfneed
    if category:
        info("limiting download to %s." % category)
        dfx=dfneed[dfneed.refseq_category == category].reset_index(drop=True)
        dffinal=pd.concat([dffinal,dfx]).drop_duplicates().reset_index(drop=True)
    info("Parsing assembly summary.")
    info("Final download count is %s" % len(dffinal))
    for a,b in dffinal.iterrows():
        url=b.ftp_path
        # name=b.organism_name
        try:
            process_url_file(url)
            TAX=b.taxid
            filename=url.split("/")[-1]+"_genomic.fna"
            get_fasta_in_kraken_format(filename, TAX,'%s_genomes.fa' % taxonomic_category)
            info("Process complete")
        except:
            info("Process failed")
    return
    

# ------------------------------------------------------------------------------
#                    MAIN FUNCTION
# ------------------------------------------------------------------------------
            
if __name__ == '__main__':
    info('Start execution', init_new_line=False)
    t0 = time.time()
    args = read_params()
    check_params(args)
    print(args)
    download_genomes_continue(args.path_to_krakendb,
                              args.taxonomic_category,
                              args.assembly_level,
                              category=args.refseq_category)

    exec_time = time.time() - t0
    info('Finish execution ({} seconds)\n'.format(round(exec_time, 2)))
    




