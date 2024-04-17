#!/bin/bash

## Get GTDB in kraken format (https://github.com/shenwei356/gtdb-taxdump)
## Scripts adopted from KMCP tutorial (https://bioinf.shenwei.me/kmcp/database/#building-databases)
## Formulated for OSCER
## Last updated Feb 22, 2024 by SJ

K2DIR=/PATH/TO/KRAKENDB

## Get NCBI Taxdump files
TAXDUMP=/PATH/TO/taxonomy
mkdir -p $TAXDUMP
cd $TAXDUMP
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvf taxdump.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gunzip nucl_gb.accession2taxid.gz

### Download and format GTDB
cd $K2DIR
wget https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz
wget https://data.gtdb.ecogenomic.org/releases/release214/214.1/bac120_metadata_r214.tsv.gz
wget https://data.gtdb.ecogenomic.org/releases/release214/214.1/ar53_metadata_r214.tsv.gz

tar -zxvf gtdb_genomes_reps.tar.gz

# rename
## Removes _genomic in string
## If you do this, change get_gtdb_in_karken_format.py to split('.fna') instead of split('_genomic')
/home/sjj9602/Programs/brename -R -p '^(\w{3}_\d{9}\.\d+).+' -r '$1.fna.gz' gtdb_genomes_reps_r214

# assembly accession -> full head
find gtdb_genomes_reps_r214/ -name "*.fna.gz" \
    | ~/Programs/rush -k 'echo -ne "{%@(.+).fna}\t$(seqkit sort --quiet -lr {} | head -n 1 | seqkit seq -n)\n" ' \
    > gtdb.name.map

cat *_metadata*.tsv \
awk -F"\t" '{OFS=FS}{print $1,$17,'
    | csvtk cut -t -f accession,gtdb_taxonomy,ncbi_taxid \
    | csvtk replace -t -p '^.._' \
    | csvtk grep -t -P <(cut -f 1 name.map) \
    | csvtk del-header \
    > taxid.map

# assembly accession -> taxid
cat accession.spec.taxid.map \
	|sed '1iaccession\tspecies\ttaxid' \
	|csvtk cut -t -f accession,taxid \
	|csvtk grep -t -P <(cut -f 1 gtdb.name.map) \
	|csvtk del-header \
	> gtdb.taxid.map

### Format GTDB for kraken2
cp gtdb.taxid.map /media/projects/k2db-GTDB/gtdb_genomes_reps_r214/database
## conda activate mpa
python ~/get_gtdb_in_kraken_format.py /media/projects/k2db-GTDB/gtdb_genomes_reps_r214/database

#######################################

## Download NCBI genomes
mkdir genomes/
cd genomes
mv ../gtdb.genomes.fa .
python update_kraken_continue.py -p $K2DIR -c fungi -a "Complete Genome" -r "representative genome"
python update_kraken_continue.py -p $K2DIR -c protozoa -a "Complete Genome" -r "representative genome"
python update_kraken_continue.py -p $K2DIR -c viral -a "Complete Genome" -r "representative genome"
# Download human genome
# Convert to kraken format
python get_fasta_in_kraken_format.py /PATH/TO/human_genome.fa 9606

### Get NCBI taxid info
## combine only fungi, viruses, protozoa, and organelles for taxid.map and name.map
## Bacteria and archaea will come from GTDB
echo -e "viral\nprotozoa\nfungi"| while read in
do 
    grep ">" $in/$in\_genomes.fa | sed 's/>//' | awk '{print $1}' | awk -F"|" '{OFS="\t"}{print $1,$3}' > $in.id.taxid.map
done
grep ">" human_genome.tax.fa |sed 's/>//g' | awk '{print $1}' | awk -F"|" '{OFS="\t"}{print $1,$3}' > human.id.taxid.map
grep ">" gtdb.genomes.fa |sed 's/>//g' | awk '{print $1}' | awk -F"|" '{OFS="\t"}{print $1,$3}' > gtdb.id.taxid.map

## Download organelle genomes
mkdir organelle
cd organelle
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.1.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.1.2.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.2.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.3.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz
gunzip *gz

grep ">" *genomic.fna |seqkit seq -n > organelle.name.map

### Get organelle taxids
awk '{print $1}' organelle.name.map > accnums.txt
## Find nucl_gb.accession2taxid in NCBI taxdump directory
## Reduce huge file size by only keeping accession numbers with underscore (refseq has underscore)
grep "_" $TAXDUMP/nucl_gb.accession2taxid > $TAXDUMP/condensed.accession2taxid

# Get accession2taxid map for my samples
python get_my_accession2taxid.py accnums.txt $TAXDUMP/condensed.accession2taxid
# output will be my.accession2taxid.map

# Get taxid map
cat my.accession2taxid.map |sed '1d' |awk -F"\t" '{OFS=FS}{print $2,$3}' >> organelle.id.taxid.map 

## Format organelle genomes for kraken format
python ~/get_organelle_in_kraken_format.py .

## Optional - For accession numbers not mapping
## Get accession numbers in 'need.txt' file
# cat need.txt |while read acc; do efetch -db nucleotide -id "$acc" -format docsum < /dev/null|xtract -pattern DocumentSummary -element AccessionVersion,TaxId >> organelle.taxid.map; done

mv organelle.id.taxid.map ../
cd ../
cat *id.taxid.map > id.taxid.map


#############################
## Create custom taxdump

## GTDB taxdump location
GTDBDUMP=gtdb-taxdump/R214.1
mkdir -p $GTDBDUMP
cd $GTDBDUMP
wget https://data.gtdb.ecogenomic.org/releases/release214/214.1/ar53_taxonomy_r214.tsv
wget https://data.gtdb.ecogenomic.org/releases/release214/214.1/bac120_taxonomy_r214.tsv
cd ..
create-taxdump $GTDBDUMP/*.tsv* --gtdb --out-dir $GTDBDUMP --force

## Get tab-separated lineage information for taxids at species level
taxonkit list --data-dir $GTDBDUMP/ --ids 1 --indent "" \
    | taxonkit filter --data-dir $GTDBDUMP/ --equal-to species \
    | taxonkit reformat --data-dir $GTDBDUMP/ --taxid-field 1 \
        --format "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" \
        -o gtdb.tsv

## NCBI taxdump tsv
# TAXDUMP refers to NCBI taxdump and is defined earlier
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
(awk '{print $_"\t"}' gtdb.tsv; cat ncbi-viral.tsv ncbi-eukaryota.tsv) \
    | taxonkit create-taxdump \
        --field-accession 1 \
        -R "superkingdom,phylum,class,order,family,genus,species,strain" \
        -O taxdump-gtdb+ncbi-eukaryote

### NEED TO CHANGE ALL TAXIDS IN NCBI GENOMES FILES TO MATCH NEW TAXDUMP
### GTDB is fine
ls *genomes.fa|grep -v "gtdb"|while read infile; do ~/fix_krakenid.py $infile ../taxdump-gtdb+ncbi-eukaryote/taxid.map; done 1> fix.o 2> fix.e &

###############################
## Get genomes used to build kraken database - build bt2 database - align

cat top_20_euk.txt |sed 's/_/ /g'|while read in; do awk -F"\t" -v a="$in" '$8 == a' ncbi-eukaryota.tsv|head -1; done > top_20_euk.ncbi.txt
cd genomes/
ls *formatted.fa|awk -F"_" '{print $1}'|while read gen; do grep ">" $gen\_genomes.formatted.fa > $gen.headers.txt; done
cat fungi.headers.txt organelle.headers.txt protozoa.headers.txt > eukaryota.headers.txt
cat fungi_genomes.formatted.fa organelle_genomes.formatted.fa protozoa_genomes.formatted.fa > eukaryota.genomes.formatted.fa
cd ../
cat top_20_euk.txt|while read in; do echo $in|sed 's/_/ /g'|grep -wf - genomes/eukaryota.headers.txt > $in.headers.txt; done
sed -i 's/>//g' *headers.txt
cat top_20_euk.txt|while read name; do seqkit grep -n -f $name.headers.txt -o $name.formatted.fa genomes/eukaryota.genomes.formatted.fa; done
