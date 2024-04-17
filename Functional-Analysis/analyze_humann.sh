#!/bin/bash

## Analyze HUMANn3.6 output

## Perform the same analyses with comparative datasets
# Those used here include Honap et al. 2023 and Velsko et al. 2017
# Move output files to respective directoryes (Genes, Pathabundance, Pathcoverage)
# Initialize variable for study name
file=CHS-Vel-WAMP

### 	Functional Analysis - including comparative datasets 	###
# Join gene family tables
# Genes is a directory - must not have anythingn but .genefamilies.tsv files here
humann_join_tables -i Genes -o $file.genefamilies.tsv
# Remove suffix from title
sed -i 's/_Abundance-RPKs//g' $file.genefamilies.tsv
# Renormalize table to copies per milion from reads per kilobase
# Normalizes for differences in depth of sequencing
humann_renorm_table -i $file.genefamilies.tsv -o $file.genefamilies.cpm.tsv -u cpm -s n
# Regroup table by KEGG orthology group using utility scripts provided by HUMANn download
# Ignore ungrouped values (-u N)
humann_regroup_table -i $file.genefamilies.cpm.tsv -o $file.ko.tsv -u N -g uniref50_ko
# Rename table using annotation mapping file provided by HUMANn
humann_rename_table -i $file.ko.tsv -n kegg-orthology -o $file.ko.named.tsv
# Split tables into stratified and unstratified
humann_split_stratified_table -i $file.ko.named.tsv -o $file.ko.named-split
sed -i 's/# Gene Family/KOname/g' $file.ko.named-split/*tsv

### Looking for herbivores and carnivores
## Muegge et al. 2011 has a list of enzyme commission numbers pertaining to diet
## used these numbers to identify KO numbers in humann output
## First get a list of just KO names (with EC)
awk -F"\t" '{print $1}' $file.ko.named_unstratified.tsv|sed '1d' > $file.kolist.txt
## Extract EC numbers and search list of KO numbers
## Output is list of KO numbers in merged table related to meat consumption
awk '{print $1}' MeatConsumption_ECnums.txt|sed '1d' |sort|uniq|while read val; do grep -w $val $file.kolist.txt; done|awk -F":" '{print $1}'|sort|uniq > $file.HerbCarn-KO.txt
## Filter merged table based on meat consumption KO values
head -1 $file.ko.named_unstratified.tsv > $file.HerbCarn.tsv; grep -f $file.HerbCarn-KO.txt $file.ko.named_unstratified.tsv >> $file.HerbCarn.tsv
## Use $file.HerbCarn.tsv and metadata file as input to HerbCarnBoxplotScript.R
## Create metadata file with SampleID in col1 and Population in col2
# data=$file.HerbCarn.tsv
# metadata=$file.metadata.txt
# featuremap=MeatConsumption_ECnums.txt

