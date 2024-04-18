# Anson-St-Ancestors-Oral-Microbiome
General workflow of bioinformatics analysis for Ansont St. Ancestor Oral Microbiome project.
Analyses were run on a combination of personal compuber (Mac), workstation (Linux), and the HPC cluster at the University of Oklahoma, OSCER.

# Table of contents

<!--ts-->
   * [Pre-processing of raw sequencing data](#pre-processing-of-raw-sequencing-data)
   * [Taxonomic classification](#taxonomic-classification)
   	  * [Download MetaPlAn4 CHOCOPhlAn database](#download-metaphlan4-chocophlan-database)
      * [Classify trimmed and merged sequences](#classify-trimmed-and-merged-sequences)
      * [Perform PERMANOVA and Biplot analysis](#perform-permanova-and-biplot-analysis)
   * [Analysis of metagenomic data in R](#analysis-of-metagenomic-data-in-r)
     * [Preparation of Phyloseq datasets](#preparation-of-phyloseq-datasets)
     * [Non-metric Multidimensional Scaling](#non-metric-multidimensional-scaling)
     * [Differential taxonomic abundances with DESeq2](#differential-taxonomic-abundances-with-deseq2)
     * [PERMANOVA and dispersion test](#permanova-and-dispersion-test)
   * [Antimicrobial resistance analysis](#antimicrobial-resistance-analysis)
     * [Preparation of the databases and Blast analysis](#preparation-of-the-databases-and-blast-analysis)
     * [AMR data analysis in R](#amr-data-analysis-in-r)    
<!--te-->

## Pre-processing of raw sequencing data
Sequencing data from the Ancestors and two negative controls were processed with AdapterRemoval to merge paired-end reads, trim Ns and poor quality reads, and eliminate fragments shorter than 30 nt, and remove adapters:

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
	--basename $sample --gzip --settings $sample.AdapterRemoval.txt
	## Combine collapsed and collapsed truncated into Analysis Ready (AR) files
	cat $sample.collapsed* > $sample.AR.fastq.gz
done
```

## Taxonomic classification 

### Download MetaPlAn4 CHOCOPhlAn database
Download a marker gene database that includes species level genome bins (SGBs) to increase classification of poorly characterized genomes. 

```bash
MPAdb=/path/to/metaphlan_databases
metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db $MPAdb
```

### Classify trimmed and merged sequences
Classify reads from Ancestors and negative controls with a minimum read length of 30 and includes unclassified sequences. Output is a profile (SAMPLE.mpa) and the bowtie2 output (SAMPLE.bt2.bz2). 

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
Perform multivariate statistical tests for metagenomes from the Ancestors:
[PCA_PERMANOVA_BIPLOT.R](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/tree/main/Taxonomic-Classification-and-Analysis/PCA_PERMANOVA_BIPLOT.R)

Perform multivariate statistical tests for metagenomes from the Ancestors and comparative populations (modern Spain, historic Asia, historic Africa, historic United Kingdom, ancient Africa, and pre-contact North America)
[PCA_PERMANOVA_BIPLOT_Comparative.R](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/tree/main/Taxonomic-Classification-and-Analysis/PCA_PERMANOVA_BIPLOT_Comparative.R)



