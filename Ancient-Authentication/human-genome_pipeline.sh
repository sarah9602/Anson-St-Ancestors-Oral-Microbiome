#!/bin/bash
## Pipeline to align analysis-ready reads to human genome.
## Script created by Sarah Johnson (2024) inspired by those from Tanvi Honap

## Install mapDamage2 (https://ginolhac.github.io/mapDamage/)

## Human genome used in Anson-Street-Ancestors project = hg19
ref=/PATH/TO/References/
cd $ref
outdir="hg19"
## Download reference and create bowtie2 index
bowtie2-build $outdir.fasta $outdir

## Assign your software variables
qualimap=/PATH/TO/qualimap_v2.2.1/qualimap
dedup=/PATH/TO/DeDup-0.12.8.jar

## Assign your personal sutdy variables
## These will change per study
study=Charleston
path=/PATH/TO/WORKING/DIRECTORY
sampdir=/PATH/TO/AnalysisReady

## Make a directory to store stats and output fastq files
cd $path
mkdir -p $outdir
cd $outdir

## Perform alignment
samp=SAMPLEID
file=$sampdir/$samp.AR.fastq.gz

## Run bowtie2 against reference
## Keep both aligned and unaligned reads to separate later.
## --no-unal flag would be used if you wanted to only keep those aligning to reference
  bowtie2 -p 5 -x $ref/$outdir -U $file -S $samp.bt2_mapped.sam
## Convert to BAM format and sort
  samtools view -bSh -@ 5 $samp.bt2_mapped.sam > $samp.bt2_mapped.bam
  samtools sort -o $samp.bt2_mapped_sorted.bam $samp.bt2_mapped.bam 
## Filter out unaligned reads
## Output is bam file with only reads aligning to human reference
## This will be used in downstream analyses.
## Suffix is reference ID
## -F is keep all except this flag - 4 is keep all but unaligned
  samtools view -b -F 4 $samp.bt2_mapped_sorted.bam > $samp.bt2_mapped_sorted_$outdir.bam
## Keep only unaligned reads
## Output will be filtered bam file with all non-human (aka microbial) reads
## This will be used in the rest of this pipeline
## -f is keep this flag - 4 is unaligned flag
  samtools view -b -f 4 $samp.bt2_mapped_sorted.bam > $samp.bt2_mapped_sorted_filtered.bam
## Convert filtered bam file to fastq
## Output will be SAMPLE.filtered.fastq.gz
  bam bam2FastQ --in $samp.bt2_mapped_sorted_filtered.bam --gzip --outBase $samp.filtered

###########
# $samp.filtered.fastq.gz used for downstream analysis
###########

## Deduplicate human read aligned bam file
  java -jar $dedup -i $samp.bt2_mapped_sorted_$outdir.bam -m -o $path/$outdir

## Remainder of pipeline only uses human read aligned files
  no_unique_reads=`samtools view -c $samp.bt2_mapped_sorted_$outdir\_rmdup.bam`

  if [ $no_unique_reads -eq 0 ]
  then
   echo "There were no unique mapped reads"

  else
  ## Run MapDamage
    mapDamage -i $samp.bt2_mapped_sorted_$outdir\_rmdup.bam -r $ref/$outdir.fasta -d $samp\_MapDamage_results

  ## Get coverage statistics
    $qualimap bamqc -bam $samp.bt2_mapped_sorted_$outdir\_rmdup.bam -outdir $samp\_QualimapStats -outformat pdf > /dev/null
  fi


## Can delete unwanted files
  # rm $samp.bt2_mapped_sorted.bam
  # rm $samp.bt2_mapped.sam
  # rm $samp.bt2_mapped.bam
  # rm $samp.filtered_1.fastq.gz
  # rm $samp.filtered_2.fastq.gz


