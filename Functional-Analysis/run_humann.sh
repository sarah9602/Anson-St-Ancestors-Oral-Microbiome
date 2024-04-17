#!/bin/bash

## Download and run HUMANn3.6 on filtered analysis-ready reads

## Classified using HUMANn v3.6
## Download HUMANn reference databases
humann_databases --download chocophlan full /PATH/TO/Humann_DB
humann_databases --download uniref uniref50_diamond /PATH/TO/Humann_DB
humann_databases --download utility_mapping full /PATH/TO/Humann_DB

## Classification takes time and memory (30-60 GB; 8-24 hr depending on metagenome size)
# Run on cluster for efficiency
## Ensure minimum read length of 30 for metaphlan functionality
# Can keep temporary folders, but accumulates much data
sample=SAMPLEID
infile=/PATH/TO/$sample.filtered.fastq.gz
WORKING=/PATH/TO/WORKING/DIRECTORY
humann -i $infile --threads 8 --search-mode uniref50 --metaphlan /PATH/TO/metaphlan/bin --metaphlan-options="--min_alignment_len 30 --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db /PATH/TO/metaphlan_databases/" --nucleotide-database /PATH/TO/Humann_DB/chocophlan --protein-database /PATH/TO/Humann_DB/uniref --output-basename $sample -o $WORKING 
# Move output files to output directories
mkdir -p Genes
mv $sample\_genefamilies.tsv Genes/
mkdir -p Pathabundance
mv $sample\_pathabundance.tsv Pathabundance/
mkdir -p Pathcoverage
mv $sample\_pathcoverage.tsv Pathcoverage/

## Optionally remove temp directory
# rm -r $WORKING/$sample\_humann_temp