# Anson-St-Ancestors-Oral-Microbiome
General workflow of bioinformatics analysis for Ansont St. Ancestor Oral Microbiome project.
Analyses were run on a combination of personal compuber (Mac), workstation (Linux), and the HPC cluster, OSCER, at the University of Oklahoma.

# Table of contents

<!--ts-->
   * [Pre-processing of raw sequencing data](#pre-processing-of-raw-sequencing-data)
   * [Taxonomic classification](#taxonomic-classification)
   	  * [Download MetaPlAn4 CHOCOPhlAn database](#download-metaplan4-chocophlan-database)
      * [Classify trimmed and merged sequences](#classify-trimmed-and-merged-sequences)
      * [Perform PERMANOVA and Biplot analysis](#perform-permanova-and-biplot-analysis)
   * [Ancient Authentication](#ancient-authentication)
     * [Perform SourceTracker2](#perform-sourcetracker2)
     * [Align trimmed and merged reads to human reference genome](#align-trimmed-and-merged-reads-to-human-reference-genome)
   * [Functional Classification](#functional-classification)
     * [Download HUMANn databases](#download-humann-databases)
     * [Classify functional potential](#classify-functional-potential)   
     * [Analyze functional profiles](#analyze-functional-profiles) 
   * [Eukaryotic Classification with Kraken2](#eukaryotic-classification-with-kraken2)
     * [Build custom database](#build-custom-database)
     * [Run Kraken and Bracken](#run-kraken-and-bracken)
     * [Analyze bracken for Eukaryotic signatures](#analyze-bracken-for-eukaryotic-signatures)

<!--te-->

## Pre-processing of raw sequencing data
Sequencing data from the Ancestors and two negative controls were processed with AdapterRemoval to merge paired-end reads, trim Ns and poor quality reads, eliminate fragments shorter than 30 nt, and remove adapters.

```bash
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
	--basename $sample --gzip --adapter-list adapters.txt --settings $sample.AdapterRemoval.txt
	## Combine collapsed and collapsed truncated into Analysis Ready (AR) files
	cat $sample.collapsed* > $sample.AR.fastq.gz
done
```

## Taxonomic classification 

### Download MetaPlAn4 CHOCOPhlAn database
Downloaded a marker gene database that includes species level genome bins (SGBs) to increase classification of poorly characterized genomes. 
Requires:
* [MetaPhlAn4](https://github.com/biobakery/MetaPhlAn)

```bash
MPAdb=/path/to/metaphlan_databases
metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db $MPAdb
```

### Classify trimmed and merged sequences
Classified reads from Ancestors and negative controls with a minimum read length of 30 and included unclassified sequences. Output was a profile (SAMPLE.mpa) and the bowtie2 alignment (SAMPLE.bt2.bz2). 

```bash
ARdir=/path/to/AnalysisReadyReads
MPAdb=/path/to/metaphlan_databases
cat $samplefile | while read sample; 
do
	metaphlan $ARdir/$sample.AR.fastq.gz --bowtie2out $sample.bt2.bz2 \
	--bowtie2db $MPAdb/ --nproc 5 --input_type fastq --read_min_len 30 -o $sample.mpa --unclassified_estimation
done
```

### Perform PERMANOVA and Biplot analysis
Performed multivariate statistical tests for metagenomes from the Ancestors:
* [PCA_PERMANOVA_BIPLOT.R](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/tree/main/Taxonomic-Classification-and-Analysis/PCA_PERMANOVA_BIPLOT.R)

Performed multivariate statistical tests for metagenomes from the Ancestors and comparative populations (modern Spain, historic Asia, historic Africa, historic United Kingdom, ancient Africa, and pre-contact North America)
* [PCA_PERMANOVA_BIPLOT_Comparative.R](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/tree/main/Taxonomic-Classification-and-Analysis/PCA_PERMANOVA_BIPLOT_Comparative.R)

## Ancient Authentication

### Perform SourceTracker2
Obtain reference shotgun metagenomes from known source environments. Ensure number of metagenomes is about the same accross all sources. Sources for this analysis - Velsko et al. 2019; HMP (Gebers et al. 2012); Sankaranarayanan et al. 2015; Rampelli et al. 2015; Obregon-Tito et al. 2015; Oh et al. 2016; Johnston et al. 2016. Create a mapping file for sources and samples (sinks) with at least 3 columns titled #SampleID, SourceSink, Env. Source samples are labeled "source" and study samples labeled "sink". Env indicates source environment, i.e. skin, ModernDentalCalculus, etc. Process source metagenomes and classify using MetaPhlAn4 with same parameters as samples. 

Update path, study, and sourcetracker directory information in the sourcetracker2 wrapper script and run.
* [sourcetracker2_metaphlan.sh](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Ancient-Authentication/sourcetracker2_metaphlan.sh)

Accompanying scripts:
* [merge_metaphlan_tables.py](https://github.com/biobakery/MetaPhlAn/blob/master/metaphlan/utils/merge_metaphlan_tables.py) (from biobakery)
* [get_sourcetracker_metaphlan_counts.py](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Ancient-Authentication/get_sourcetracker_metaphlan_counts.py)
* [merge_metaphlan_sourcetracker_counts.py](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Ancient-Authentication/merge_metaphlan_sourcetracker_counts.py)
* [sourcetracker_feature_contributions.py ](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Ancient-Authentication/sourcetracker_feature_contributions.py) (Optional)

Visualize SourceTracker2 output as pie plots:
* [Rscript_PiePlots_SourceTracker2.R](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Ancient-Authentication/Rscript_PiePlots_SourceTracker2.R)

### Align trimmed and merged reads to human reference genome
Reads were aligned to the human genome (hg19) to assess endogenous human DNA content and filter human DNA from the Ancestors' metagenomes. Resulting filtered.fastq.gz files were used as input for other downstream analyses.
Reference for current project was [hg19](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.13/).

Requires:
* [MapDamage2](https://ginolhac.github.io/mapDamage/)

Update the following and run.
* [human-genome_pipeline.sh](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Ancient-Authentication/human-genome_pipeline.sh)


## Functional Classification

### Download HUMANn databases
Downloaded full CHOCOPhlAn and UniRef50 databases along with utility scripts.
Reuqires at least [HUMANn v3.6](https://github.com/biobakery/humann)

```bash
humann_databases --download chocophlan full /PATH/TO/Humann_DB
humann_databases --download uniref uniref50_diamond /PATH/TO/Humann_DB
humann_databases --download utility_mapping full /PATH/TO/Humann_DB
```

### Classify functional potential
Classified filtered metagenomes using HUMANn3.6. Ensure minimum read length of 30 for metaphlan functionality. Classification takes time and memory (30-60 GB; 8-24 hr depending on metagenome size). Ran on a HPC cluster for efficiency. You can keep temporary folders, but this accumulates much data.
 ```bash
sample=SAMPLEID
infile=/PATH/TO/$sample.AR.fastq.gz
WORKING=/PATH/TO/WORKING/DIRECTORY
humann -i $infile --threads 8 --search-mode uniref50 --metaphlan /PATH/TO/metaphlan/bin --metaphlan-options="--min_alignment_len 30 --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db /PATH/TO/metaphlan_databases/" --nucleotide-database /PATH/TO/Humann_DB/chocophlan --protein-database /PATH/TO/Humann_DB/uniref --output-basename $SAMPLE -o $WORKING 

# Move output files to output directories
mkdir Genes
mv $sample\_genefamilies.tsv Genes/
mkdir Pathabundance
mv $sample\_pathabundance.tsv Pathabundance/
mkdir Pathcoverage
mv $sample\_pathcoverage.tsv Pathcoverage/

## Optionally remove temp directory
# rm -r $WORKING/$sample\_humann_temp
 ```

 ### Analyze functional profiles
Perform the same analyses with comparative datasets. Those used here include Honap et al. 2023 and Velsko et al. 2017.

```bash
# Initialize variable for study name
file=CHS-Vel-WAMP

# Join gene family tables
humann_join_tables -i Genes -o $file.genefamilies.tsv
# Remove suffix from title
sed -i 's/_Abundance-RPKs//g' $file.genefamilies.tsv
# Renormalize table to copies per milion from reads per kilobase
# Normalize for depth of sequencing differences
humann_renorm_table -i $file.genefamilies.tsv -o $file.genefamilies.cpm.tsv -u cpm -s n
# Regroup table by KEGG orthology group using utility scripts provided by HUMANn download
# Ignore ungrouped values (-u N)
humann_regroup_table -i $file.genefamilies.cpm.tsv -o $file.ko.tsv -u N -g uniref50_ko
# Rename table using annotation mapping file provided by HUMANn
humann_rename_table -i $file.ko.tsv -n kegg-orthology -o $file.ko.named.tsv
# Split tables into stratified and unstratified
humann_split_stratified_table -i $file.ko.named.tsv -o $file.ko.named-split
sed -i 's/# Gene Family/KOname/g' $file.ko.named-split/*tsv
```
Meat consumption analysis performed using a list of enzyme commision numbers outlined by Muegge et al. (2011). These numbers were used to identify KEGG orthologous numbers in the HUMANn output.

Required file: 
* [MeatConsumption_ECnums.txt](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Files/MeatConsumption_ECnums.txt)

```bash
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
```

Visualize comparison with the following script:
* [HerbCarnBoxplotScript.R](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Functional-Analysis/HerbCarnBoxplotScript.R)


## Eukaryotic Classification with Kraken2

### Build custom database
Database was constructed using GTDB bacteria and archaea, NCBI viruses, fungi, protozoa, plastids, and chloroplasts.
Steps for constructing the custom library and building both the kraken and bracken databases are outlined in the following script:
* [custom_krakendb.sh](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Kraken-Analysis/custom_krakendb.sh)

Requires:
* [Kraken2](https://github.com/DerrickWood/kraken2)
* [Bracken](https://github.com/jenniferlu717/Bracken)
* [taxonkit](https://bioinf.shenwei.me/taxonkit/)
* [brename](https://github.com/shenwei356/brename)
* [rush](https://github.com/shenwei356/rush)
* [Biopython](https://biopython.org/wiki/Download)

Accompanying scripts:
* [get_gtdb_in_kraken_format.py](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Kraken-Analysis/get_gtdb_in_kraken_format.py)
* [update_kraken_continue.py](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Kraken-Analysis/update_kraken_continue.py)
* [get_my_accession2taxid.py](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Kraken-Analysis/get_my_accession2taxid.py)
* [get_organelle_in_kraken_format.py](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Kraken-Analysis/get_organelle_in_kraken_format.py)
* [fix_krakenid.py](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Kraken-Analysis/fix_krakenid.py)
* [get_fasta_in_kraken_format.py](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Kraken-Analysis/get_fasta_in_kraken_format.py)

### Run Kraken and Bracken
Classify Ancestors' metagenomes using custom Kraken2 database and re-estiimate abundances using Bracken for the purposes of identifying reads attributed to eukaryotic sources, or dietary constituents. Running Kraken requires ~ 450GB of memory.

```bash
K2DB=/path/to/k2db_GTDB
samplefile=/path/to/file/with/sample/ids
FILTdir=/path/to/FilteredAnalysisReady_reads
reportdir=/PATH/TO/OUTPUT/REPORTS

## Create directory for taxonomic levels to summarize with bracken
mkdir -p $reportdir/bracken-G
mkdir -p $reportdir/bracken-S

## Database with GTDB bacteria and archaea, NCBI viruses, fungi, protozoa, plastids, and chloroplasts
## Built and run with default parameters
## Requires ~ 450 GB memory
cat $samplefile|while read sample;
do
	## Run Kraken2
	kraken2 --db $K2DB --threads 20 --report $reportdir/$sample.k2report --report-minimizer-data --gzip-compressed --classified-out $reportdir/$sample.classified.fastq --output $reportdir/$sample.k2out $FILTdir/$sample.filtered.fastq.gz
	## If running with --report-minimizer-data, these two columns should be removed before running bracken
	## remove minimizer data
	cat $reportdir/$sample.k2report|cut -f4,5 --complement > $reportdir/$sample.k2
	## Run Bracken
	bracken -d $K2DB -i $reportdir/$sample.k2 -o $reportdir/bracken-S/$sample.S.bracken -w $reportdir/bracken-S/$sample.S.k2report.bracken -r 35 -l S -t 10
	bracken -d $K2DB -i $reportdir/$sample.k2 -o $reportdir/bracken-G/$sample.G.bracken -w $reportdir/bracken-G/$sample.G.k2report.bracken -r 35 -l G -t 10
done
```

### Analyze bracken for Eukaryotic signatures
Determine if DNA from eukaryotes were identified in metagenomic samples as possible direct evidence of dietary constituents. To assess validity of diety through this avenue of exploration, reads must be mapped to a reference and damage patterns must be assessed. 

Requires:
* [KrakenTools](https://github.com/jenniferlu717/KrakenTools)
* [MapDamage2](https://ginolhac.github.io/mapDamage/)

```bash
K2DB=/path/to/k2db_GTDB
reportdir=/PATH/TO/OUTPUT/REPORTS
cd $reportdir

## Convert desired taxonomic level to metaphlan style
ls *S.k2report.bracken|while read in; do 
	samp=`echo $in|awk -F"." '{print $1}'`; \
	kreport2mpa.py -r $in -o k2report-mpa/$samp.k2.mpa --display-header; \
done

combine_mpa.py -i *mpa -o Charleston_CombinedKraken2mpa.s.txt

## Get eukaryotes
head -1 Charleston_CombinedKraken2mpa.s.txt > Charleston_Eukaryote.s.txt
grep "Eukaryo" Charleston_CombinedKraken2mpa.s.txt |grep "s__" >> Charleston_Eukaryote.s.txt
## Rremove suffix from header
sed -i 's/.S.k2report.bracken//g' Charleston_Eukaryote.s.txt

## Sort and filter combined Eukaryotic table and get top 20 eukaryotes
python get_kraken_eukaryotes.py Charleston_Eukaryote.s.txt
## output will be Charleston_Eukaryote.s.top20.txt and top_20_euk.txt file containing list of top taxa

cd $K2DB/genomes
cat fungi_genomes.formatted.fa organelle_genomes.formatted.fa protozoa_genomes.formatted.fa > eukaryota.genomes.formatted.fa
grep ">" eukaryota.genomes.formatted.fa|sed 's/>//g' > eukaryota.headers.txt
cd $reportdir
## Get top20 eukaryotes header file
cat top_20_euk.txt|while read in; do \
	echo $in|sed 's/_/ /g'|grep -wf - $K2DB/genomes/eukaryota.headers.txt > $in.headers.txt; \
done
## Use header file to extract reads from formatted genome files
cat top_20_euk.txt|while read name; do \
	seqkit grep -n -f $name.headers.txt -o $name.formatted.fa $K2DB/genomes/eukaryota.genomes.formatted.fa; \
done
```

Align sample reads to new genome files using the following:
* [kraken_euk_mapping.sh](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Kraken-Analysis/kraken_euk_mapping.sh).



