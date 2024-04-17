#!/bin/bash

## Download databases and run MetaPhlAn4

## Download desired MetaPhlAn4 database or let download occur automatically
MPAdb=/path/to/metaphlan_databases
metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db $MPAdb

## Classified using MetaPhlAn v4.0.1
## Ensure minimum read length of 30
## Ensure unclassified_estimation to obtain overall relative abundance

# Set variables
samplefile=/path/to/file/with/sample/ids
ARdir=/path/to/AnalysisReadyReads

cat $samplefile | while read sample; 
do
	metaphlan $ARdir/$sample.AR.fastq.gz --bowtie2out $sample.bt2.bz2 \
	--bowtie2db $MPAdb/ --nproc 5 --input_type fastq --read_min_len 30 -o $sample.mpa --unclassified_estimation
done

# Output will be SAMPLE.mpa and SAMPLE.bt2.bz2