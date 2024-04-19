#!/bin/bash
#
############

## Script for running SourceTracker2 from metaphlan output
# Ensure sources and samples (sinks) classified with the same database and metaphlan version using same parameters.
# Blanks can be used as sources if species are identified.
# Rarefying to 20K is arbitrary, but the higher you go, the longer it takes.
# Resulting mixing_proportions.txt file can be used as input for R to visualize plots.
# Metaphlan mapping file noted as otu-map-TAXONLEVEL.txt are arbitrarily created below. 
# Use map file specific for taxon level.

## Created by Sarah Johnson (2024)

## Need sourcetracker2 environment
conda create -n st2 -c conda-forge python=3.8 numpy scipy scikit-bio h5py hdf5 seaborn biom-format
conda activate st2
pip install https://github.com/biota/sourcetracker2/archive/master.zip

#############
###		Set variables		###
path=/PATH/TO/METAPHLAN/OUTPUT; \
study=Charleston; \
## Idntify sourcetracker path
SoTr=/PATH/TO/SourceTracker/DATA; \
cd $path
# Taxon level can be phylum (p), class (c), order (o), family (f), genus (g), or species (s)
taxonlevel='s'
# Not tax is the level just below desired level to remove this from output
# 't' is strain level in MetaPhlAn4
nottax='t'

###############################################
###		Make Mapping File 		###
# No underscores in sample titles!
# Create a mapping file for sources and samples (sinks) with at least 3 columns titled #SampleID, SourceSink, Env
# Blanks should be labeled as "source" if species are identified
# File is ST-$study.MappingFile.txt

###############################################
### 	DO BEFORE ANALYSIS		###
# Use merge_metaphlan_tables.py provided by metaphlan
# Must run metaphlan with --unclassified_estimation flag
# Cannot have taxid column (column 2)
# sample IDS must match mapping file
# remove suffix if necessary


###############################################
###		Process SourceTracker 		###
## Create mapping file to arbitrary OTUID at desired taxonlevel
# Get species.txt file from metaphlan release
awk '{print $2}' /PATH/TO/METAPHLAN/DATABASES/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt|sed 's/,/\t/g'|while read line; do for l in $line; do echo $l; done; done|sort|uniq|awk -F"\|$nottax\_" '{print $1}'|awk -F"\t" '{OFS=FS}{print "Taxon"++i,$0}' > $SoTr/Metaphlan/otu-map-$taxonlevel.txt

## Merge Source metaphlan output files
# Metaphlan output in Metaphlan directory
merge_metaphlan_tables.py -o $SoTr/Metaphlan/merged.ST.txt $SoTr/Metaphlan/*mpa

## Get a mapping file with sample names in col1 and read counts in col2
ls $SoTr/Metaphlan/*mpa|while read in; do \
    samp=`echo $in|awk -F"/" '{print $NF}'|awk -F"." '{print $1}'`; \
    ct=`grep "reads processed" $in|sed 's/#//'|awk '{print $1}'|bc`; \
    echo -e "$samp\t$ct"; done > $SoTr/Metaphlan/ST.counts.txt

## Output is merged metaphlan table with counts calculated using total reads and relative abundances
# Taxid column omitted
python scripts/get_sourcetracker_metaphlan_counts.py $SoTr/Metaphlan/merged.ST.txt $SoTr/Metaphlan/ST.counts.txt
## Output is $SoTr/Metaphlan/merged.ST.counts.txt and will be combined to merged sample output.

###############################################
### 	Process Samples 	###
# Cannot have taxid column (column 2)
merge_metaphlan_tables.py -o merged.$study.txt *mpa
# sample IDS must match mapping file
# remove suffix if necessary

## Get count mapping file as above
ls *mpa|while read in; do samp=`echo $in|awk -F"." '{print $1}'`; \
	ct=`grep "reads processed" $in|sed 's/#//'|awk '{print $1}'|bc`; \
    echo -e "$samp\t$ct"; done > $study.counts.txt

python scripts/merge_metaphlan_sourcetracker_counts.py -s merged.$study.txt -st $SoTr/Metaphlan/merged.ST.counts.txt -m $study.counts.txt

## Summarize MetaPhlAn output to desired taxonomic level
head -1 ST-merged.$study.counts.txt > ST-merged.$study.counts.$taxonlevel.biom.txt; cat ST-merged.$study.counts.txt|grep "$taxonlevel\__"|grep -v "$nottax\__" >> ST-merged.$study.counts.$taxonlevel.biom.txt
## Add taxonomic lineage to final column in temporary file
awk -F'\t' '{OFS=FS}{print $0,$1}' ST-merged.$study.counts.$taxonlevel.biom.txt|sed 's/clade_name/taxonomy/g' > tmp.txt
## Convert first column to OTUID given a mapping file
awk -F"\t" -v OFS="\t" 'FNR == NR{a[$2] = $1; next} {$1 =a[$1]}1' $SoTr/Metaphlan/otu-map-$taxonlevel.txt tmp.txt > ST-merged.$study.$taxonlevel.counts.withtaxon.txt
## Remove temporary file.
rm tmp.txt
## Convert to biom file
biom convert -i ST-merged.$study.$taxonlevel.counts.withtaxon.txt -o ST-merged.$study.$taxonlevel.OTU.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
## Summarize taxa
biom summarize-table -i ST-merged.$study.$taxonlevel.OTU.biom -o ST-merged.$study.$taxonlevel.OTU.biom_summary.txt

##################################
### Run SourceTracker2
conda activate st2
## Rarefy using sourcetracker2 script
## Arbitrarily chose 20K - much higher values significanlty slow down the program.
rare=20000
## Run sourcetracker for desired taxonomic level
# Takes time depending on rarefaction
# Rarefy both samples and sources to same level
# Alter alpha values to as needed
# per_sink_feature_assignments will output feature_table for each sample detailing contributing OTU ids
sourcetracker2 gibbs -i ST-merged.$study.$taxonlevel.OTU.biom -m ST-$study.MappingFile.txt -o sourcetracker2.$rare.$taxonlevel --source_rarefaction_depth $rare --sink_rarefaction_depth $rare  --jobs 5 --alpha1 0.01 --alpha2 1.0 --per_sink_feature_assignments 1> st2.$taxonlevel.log 2> st2.$taxonlevel.err &

##################################

## Script to convert biom to tsv
## Useful if you want to inspect your input biom file
# biom convert -i ST-merged.$study.$taxonlevel.OTU.biom -o ST-merged.$study.$taxonlevel.OTU.biom.txt --to-tsv --header-key taxonomy --output-metadata-id "#OTUID"

##################################

## Inspect top contributing taxa (optional)
sourcetracker_feature_contributions.py -s /PATH/TO/sourcetracker2.$rare.$taxonlevel -st $SoTr/Metaphlan -t $taxonlevel
