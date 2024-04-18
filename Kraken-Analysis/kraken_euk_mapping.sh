#!/bin/bash
## This pipeline will align samples to genomes used to feed the kraken database
## Requires mapDamage to be in path
## Performed with BWA v0.7.15-r1140
## Script written by Sarah Johnson (2024) inspired by Tanvi Honap.

## Set variables
reportdir=/PATH/TO/OUTPUT/REPORTS
path=/PATH/TO/WORKING/DIR
FILTdir=/path/to/FilteredAnalysisReady_reads
VarScan=/PATH/TO/VarScan.v2.3.9.jar
qualimap=/PATH/TO/qualimap_v2.2.1/qualimap
dedup=/PATH/TO/DeDup-0.12.8.jar

cd $reportdir
## Index these fasta files
ls *formatted.fa|whlie read in; \
do bwa index $in; \
done

## Get a list of samples to align to each one
cat top_20_euk.txt |\
while read name; \
do echo $name|sed 's/_/ /g'|grep -wf - *k2|awk -F"." '{print $1}' \
> $name.samples.txt; done

ls *samples.txt| awk -F"." '{print $1}'|\
while read name; do \
mkdir -p $name; \
echo -e "Sample\tUnique mapped reads\tAverage length of mapped reads\tMean coverage\tStd dev of mean coverage\t% reference covered >= 1X\t% reference covered >= 5X\tNumber of bp in reference" > $name.pipeline_summary.txt; \
cat $name.samples.txt| \
while read samp; do \
## reference is eukaryote name AKA name
## Query is filtered fastq.gz file
## Output will be in created directories
bwa aln -l 1024 -q 37 -n 0.1 -t 5 $name.formatted.fa $FILTdir/$samp.filtered.fastq.gz > $name/$samp.sai
bwa samse $name.formatted.fa  $name/$samp.mapped.sai $FILTdir/$samp.filtered.fastq.gz > $name/$samp.sam
samtools view -bSh -@ 5 $name/$samp.sam > $name/$samp.bam
## Filter out unmappped and low quality reads
## samtools view -bh -F4 -q displays the previous output as a BAM file (b), and
## includes the header (h), and skips alignments containing the 4 flag (0x4 segment unmapped)
samtools view -bh -F4 $name/$samp.bam > $name/$samp.mapped.bam

## Sort alignments by leftmost coordinates
samtools sort -o $name/$samp.mapped_sorted.bam $name/$samp.mapped.bam 
java -jar $dedup -i $name/$samp.mapped_sorted.bam -m -o $name
no_unique_reads=`samtools view -c $name/$samp.mapped_sorted_rmdup.bam`
if [ $no_unique_reads -eq 0 ]
  then
   length="There were no unique mapped reads"
   mean_cov="N/A"
   sd_cov="N/A"
   cov1="N/A"
   cov5="N/A"
   refbp="N/A"

  else
  	length=`samtools view $name/$samp.mapped_sorted_rmdup.bam |awk '{SUM+=length($10);DIV++}END{print SUM/DIV}'`
  	## Run mapdamage
  	mapDamage -i $name/$samp.mapped_sorted_rmdup.bam -r $name.formatted.fa -d $name/$samp.MapDamage_results 
    mv $name/$samp.MapDamage_results/Fragmisincorporation_plot.pdf $name/$name.$samp.Fragmisincorporation_plot.pdf
    ## Get quality statistics of alignment
    $qualimap bamqc -bam $name/$samp.mapped_sorted_rmdup.bam -outdir $name/$samp.QualimapStats -outformat pdf > /dev/null; \

    if [ $? -eq 0 ]; \
    then \
     mean_cov=`grep "mean coverageData" $name/$samp.QualimapStats/genome_results.txt | cut -d "=" -f2`; \
     sd_cov=`grep "std coverageData" $name/$samp.QualimapStats/genome_results.txt | cut -d "=" -f2`; \
     cov1=`grep "reference with a coverageData >= 1X" $name/$samp.QualimapStats/genome_results.txt | awk '{print $4}'`; \
     cov5=`grep "reference with a coverageData >= 5X" $name/$samp.QualimapStats/genome_results.txt | awk '{print $4}'`; \
     refbp=`grep "number of bases" $name/$samp.QualimapStats/genome_results.txt |awk '{print $5}'|sed 's/,//g'|bc`; \
   
    else \
     mean_cov="N/A"; \
     sd_cov="N/A"; \
     cov1="N/A"; \
     cov5="N/A"; \
     refbp="N/A"; \
 fi; \
fi; \

echo -e "$samp\t$no_unique_reads\t$length\t$mean_cov\t$sd_cov\t$cov1\t$cov5\t$refbp" >> $name.pipeline_summary.txt; \
done; \
done


