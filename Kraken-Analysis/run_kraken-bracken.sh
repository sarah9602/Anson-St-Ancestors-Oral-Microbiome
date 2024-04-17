#!/bin/bash

## Run kraken and bracken on samples
K2DB=/path/to/k2db_GTDB
samplefile=/path/to/file/with/sample/ids
FILTdir=/path/to/FilteredAnalysisReady_reads
reportdir=/PATH/TO/OUTPUT/REPORTS

## Create directory for taxonomic levels to summarize with bracken
mkdir -p $reportdir/bracken-G
mkdir -p $reportdir/bracken-S

## Database with GTDB bacteria and archaea, NCBI viruses, fungi, protozoa, plastids, and chloroplasts
## Built with default parameters
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