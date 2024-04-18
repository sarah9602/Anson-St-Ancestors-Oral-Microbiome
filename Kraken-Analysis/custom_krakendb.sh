#!/bin/bash

## Create custom kraken2 database
## Script written by Sarah Johnson (2024)
## Database includes GTDB bacteria and archaea, NCBI fungi, protozoa, viruses, mitochondria, and plastids
## Requirements:
##		taxonkit (https://bioinf.shenwei.me/taxonkit/)
## 		brename (https://github.com/shenwei356/brename)
##		rush (https://github.com/shenwei356/rush)
##		biopython (https://biopython.org/wiki/Download)
## Some scripts inspired by gtdb-taxdump (https://github.com/shenwei356/gtdb-taxdump)
## Some scripts adopted from KMCP tutorial (https://bioinf.shenwei.me/kmcp/database/#building-databases)

## Download NCBI taxdump files
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir -p NCBI-taxonomy
mv taxdump.tar.gz NCBI-taxonomy/
tar -xvf NCBI-taxonomy/taxdump.tar.gz

TAXDUMP=/PATH/TO/NCBI-taxonomy
K2DB=/path/to/k2db_GTDB
mkdir $K2DB
cd $K2DB

##########################
### FORMAT GTDB
## Using GTDB version R214.1
wget https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz

mkdir -p gtdb
tar -zxvf gtdb_genomes_reps.tar.gz

# rename
## Removes _genomic in string
brename -R -p '^(\w{3}_\d{9}\.\d+).+' -r '$1.fna.gz' gtdb

# assembly accession -> full head
find gtdb/ -name "*.fna.gz" \
    | rush -k 'echo -ne "{%@(.+).fna}\t$(seqkit sort --quiet -lr {} | head -n 1 | seqkit seq -n)\n" ' \
    > gtdb.name.map

# assembly accession -> taxid
cat accession.spec.taxid.map \
	|sed '1iaccession\tspecies\ttaxid' \
	|csvtk cut -t -f accession,taxid \
	|csvtk grep -t -P <(cut -f 1 gtdb.name.map) \
	|csvtk del-header \
	> gtdb.taxid.map

### Format GTDB for kraken2
cp gtdb.taxid.map $K2DB/gtdb_genomes_reps_r214/database
## conda activate mpa
python ~/get_gtdb_in_kraken_format.py $K2DB/gtdb_genomes_reps_r214/database

#######################################
### FORMAT NCBI
## combine only fungi, viruses, protozoa, and organelles for taxid.map and name.map
## Bacteria and archaea will come from GTDB

echo -e "fungi\nprotozoa\nviral"|while read in; \
do \
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/$in/assembly_summary.txt; \
mv assembly_summary.txt $in.assembly_summary.txt; \
cat $in.assembly_summary.txt|sed '1,2d'|cut -f 1,6 > $in.taxid.map; \
cat $in.assembly_summary.txt|sed '1,2d'|cut -f 1,8 > $in.name.map; \
cat $in.taxid.map | csvtk freq -Ht -f 2 -nr | taxonkit lineage -r -n -L --data-dir $TAXDUMP/ | taxonkit reformat -I 1 -f '{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}' --data-dir $TAXDUMP/ | csvtk add-header -t -n 'taxid,count,name,rank,superkindom,phylum,class,order,family,genus,species' > refseq-$in.taxid.map.stats.tsv; \
taxonkit reformat --data-dir $TAXDUMP/ --taxid-field 2 --pseudo-strain --format "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" $in.taxid.map |csvtk cut -H -t -f -2 -o refseq-$in.tsv; \
done

## Download and format NCBI genomes
## -a provide space separated values
##		if space deviding term (like Complete Genome), provide in quotes
### usage: update_kraken_continue.py [-h] [-p PATH_TO_KRAKENDB] [-c TAXONOMIC_CATEGORY] [-a ASSEMBLY_LEVEL [ASSEMBLY_LEVEL ...]] [-r REFSEQ_CATEGORY]

### optional arguments:
###   -h, --help            show this help message and exit
###   -p PATH_TO_KRAKENDB, --path_to_krakendb PATH_TO_KRAKENDB
###                         Path to where the kraken db will be built. If it doesn't exist, it will be created. Default is the working directory.
###   -c TAXONOMIC_CATEGORY, --taxonomic_category TAXONOMIC_CATEGORY
###                         Taxonomic domain to download (Required). Options include: bacteria, fungi, protozoa, archaea, viral
###   -a ASSEMBLY_LEVEL [ASSEMBLY_LEVEL ...], --assembly_level ASSEMBLY_LEVEL [ASSEMBLY_LEVEL ...]
###                         List of desired assembly completion levels (Optional). Options include: Complete Genome (default), Chromosome, Scaffold, Contig
###   -r REFSEQ_CATEGORY, --refseq_category REFSEQ_CATEGORY
###                         Level of assembly completion (Optional). Options include: representative genome, reference genome Take care when choosing.
python ~/update_kraken_continue.py -p $K2DB/genomes -c fungi -a "Complete Genome" -r "representative genome" 
python ~/update_kraken_continue.py -p $K2DB/genomes -c viral -a "Complete Genome" -r "representative genome" 
python ~/update_kraken_continue.py -p $K2DB/genomes -c protozoa -a "Complete Genome" -r "representative genome" 
mv genomes/*genomes.fa .

wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.1.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.1.2.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.2.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.2.2.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.3.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.1.genomic.fna.gz
mv *genomic.fna.gz $K2DB/genomes/
cd $K2DB/genomes
gunzip *genomic.fna.gz

echo -e "1\n2\n3"|while read num; do \
grep ">" plastid.$num.1.genomic.fna |seqkit seq -n >> plastid.name.map; \
done
echo -e "1\n2"|while read num; do \
grep ">" plastid.$num.2.genomic.fna |seqkit seq -n >> plastid.name.map; \
done

grep ">" mitochondrion.1.1.genomic.fna |seqkit seq -n > mitochondrion.name.map; \
cat mitochondrion.name.map plastid.name.map > organelle.name.map

### Get organelle taxids
awk '{print $1}' organelle.name.map > accnums.txt
## Find nucl_gb.accession2taxid in NCBI taxdump directory
## Reduce huge file size by only keeping accession numbers with underscore (refseq has underscore)
grep "_" nucl_gb.accession2taxid > condensed.accession2taxid

### Python script for mapping organelle accession numbers to mapping file
get_my_accession2taxid.py accnums.txt condensed.accession2taxid

## Format organelle genomes for kraken format
python ~/get_organelle_in_kraken_format.py $K2DB/genomes

#############################
## Create custom taxdump

## Get tab-separated lineage information for taxids at species level

## Create GTDB taxdump as shown here (https://github.com/shenwei356/gtdb-taxdump)
taxonkit list --data-dir gtdb-taxdump/R214.1/ --ids 1 --indent "" \
    | taxonkit filter --data-dir gtdb-taxdump/R214.1/ --equal-to species \
    | taxonkit reformat --data-dir gtdb-taxdump/R214.1/ --taxid-field 1 \
        --format "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" \
        -o gtdb.tsv

## Viruses taxid 10239
taxonkit list --data-dir $TAXDUMP/ --ids 10239 --indent "" \
    | taxonkit filter --data-dir $TAXDUMP/ --equal-to species --lower-than species \
    | taxonkit reformat --data-dir $TAXDUMP/ --taxid-field 1 \
        --pseudo-strain --format "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" \
        -o ncbi-viral.tsv

## Eukaryota taxid 2759
taxonkit list --data-dir $TAXDUMP/ --ids 2759 --indent "" \
    | taxonkit filter --data-dir $TAXDUMP/ --equal-to species --lower-than species \
    | taxonkit reformat --data-dir $TAXDUMP/ --taxid-field 1 \
        --pseudo-strain --format "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" \
        -o ncbi-eukaryota.tsv

## Create taxdump
## we use --field-accession  1 to output the mapping file between old taxids and new ones.
(awk '{print $_"\t"}' gtdb.tsv; cat ncbi-viral.tsv) \
    | taxonkit create-taxdump \
        --field-accession 1 \
        -R "superkingdom,phylum,class,order,family,genus,species,strain" \
        -O taxdump-gtdb+ncbi-eukaryote

## Fix genome headers to match new taxdump
## Do not need to fix gtdb genomes
cd genomes
ls *genomes.fa|grep -v "gtdb"|while read infile; do ~/fix_krakenid.py $infile /PATH/TO/taxdump-gtdb+ncbi-eukaryote/taxid.map; done

## Download human genome T2T-CHM13v2.0 name it human_genome.fa
## Format it for kraken and replace taxid to match new taxdump
## NCBI taxid for Homo sapiens is 9606
newtax=`grep -w 9606 /PATH/TO/taxdump-gtdb+ncbi-eukaryote/taxid.map|awk '{print $2}'
get_fasta_in_kraken_format.py human_genome.fa $newtax

#############################
## Add genomes to kraken library
kraken2-build --add-to-library human_genome.tax.fa --threads 4 --db $K2DB
kraken2-build --add-to-library viral_genomes.formatted.fa --threads 4 --db $K2DB
kraken2-build --add-to-library fungi_genomes.formatted.fa --threads 4 --db $K2DB
kraken2-build --add-to-library organelle_genomes.formatted.fa --threads 4 --db $K2DB
kraken2-build --add-to-library protozoa_genomes.formatted.fa --threads 4 --db $K2DB
kraken2-build --add-to-library gtdb.genomes.fa --threads 4 --db $K2DB

#############################
## Build Kraken2 DB
# Requires ~ 400 GB memory
kraken2-build --build --threads 24 --db $K2DB

# Inspect database
kraken2-inspect --db $K2DB --threads 24 > $K2DB/$K2DB.inspect.k2

#############################
## Build Bracken DB
# May require ~ 400 GB memory
bracken-build -t 20 -d $K2DB -k 35 -l 35
