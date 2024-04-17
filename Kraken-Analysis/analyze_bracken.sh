#!/bin/bash
## Analyze bracken2 output for Eukaryotic signatures

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

## Align sample reads to new genome files using kraken_euk_mapping.sh script