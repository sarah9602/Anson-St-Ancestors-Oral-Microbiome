#!/bin/bash

## Identify Adapters
## Optional - Illumina adapters used
path=/path/to/folder/with/R1-R2/files
samplefile=/path/to/file/with/sample/ids.txt

cat $samplefile | while read sample; 
do
	AdapterRemoval --identify-adapters --file1 $path/$sample.R1.fastq.gz --file2 $path/$sample.R2.fastq.gz --threads 2 1> $sample.AR.log 2> $sample.AR.er
done

# Obtain adapter list
for i in *AR.log; 
do
	adpt1=`grep adapter1 $i|awk '{print $NF}'|awk -F"NNN" '{print $1}'`
	len1=${#adpt1}; cons1=`grep -A 2 adapter1 $i|tail -1|awk '{print $NF}'`
	string1=${cons1:0:len1}; string2=`grep -A 2 adapter2 $i|tail -1|awk '{print $NF}'|awk -F"GTG" '{OFS=FS}{print $1,$2}'`
	echo -e "$string1\t$string2"
done > adapters.txt


## Adapter Removal
## Trim reads with Phred quality < 30 and a length < 30 nt
## Trim Ns
## Collapse reads with minimum alignment length of 10
cat $samplefile | while read sample; 
do 
	AdapterRemoval --threads 2 --file1 $path/$sample.R1.fastq.gz --file2 $path/$sample.R2.fastq.gz \
	--trimns --trimqualities --minalignmentlength 10 --collapse --minquality 30 --minlength 30 \
	--outputcollapsed $sample.collapsed.fastq.gz --outputcollapsedtruncated $sample.collapsed.truncated.fastq.gz \
	--basename $sample --gzip --settings $sample.AdapterRemoval.txt
	## Combine collapsed and collapsed truncated into Analysis Ready (AR) files
	cat $sample.collapsed* > $sample.AR.fastq.gz
done