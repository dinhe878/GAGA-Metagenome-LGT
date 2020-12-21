#!/bin/bash

#########################################################
# loading computerome modules
#########################################################

module load tools perl edirect/7.50 mmseqs2/release_12-113e3 emboss/6.6.0

#########################################################
# prepare taxidmap for mmseq database
#########################################################

# NCBI API KEY
export NCBI_API_KEY=593466ec18bd4cdf7c4cfd2ea54ff3908209
# location of BLASTdb
DB=/home/people/dinghe/ku_00039/people/dinghe/BLASTdb/mmseq/

# getting NCBI taxonomy
cd $DB
mkdir ncbi-taxdump && cd ncbi-taxdump
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xzvf taxdump.tar.gz
cd -

# preparing prok database
# genomes
grep ">" $DB/DB_source_fasta/prok_1908_genomes.fa | awk -F ' ' '{ gsub(">","",$1); print $1 }' > prok_1908_genomes.ids
cat prok_1908_genomes.ids | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element AccessionVersion,TaxId > prok_1908_genomes.taxidmapping
mmseqs createdb $DB/DB_source_fasta/prok_1908_genomes.fa mmseq.genome.prokDB
mmseqs createtaxdb mmseq.genome.prokDB tmp --ncbi-tax-dump ncbi-taxdump --tax-mapping-file prok_1908_genomes.taxidmapping

# proteomes (epost pipe didn't finish probably due to server issues)
mmseqs createdb $DB/DB_source_fasta/prok_1908_proteomes.fa mmseq.proteome.prokDB
mmseqs createindex mmseq.proteome.prokDB tmp
#grep ">" $DB/DB_source_fasta/prok_1908_proteomes.fa | awk -F ' ' '{ gsub(">","",$1); print $1 }' | awk -F '_' '{ print $3 }' > prok_1908_proteomes.ids
#cat prok_1908_proteomes.ids | epost -db protein | esummary -db nuccore | xtract -pattern DocumentSummary -element AccessionVersion,TaxId > prok_1908_proteomes.taxidmapping


# preparing insect database (somehow not working)
#grep ">" $DB/DB_source_fasta/insect_43_genomes.fa | awk -F ' ' '{ gsub(">","",$1); print $1 }' > insect_43_genomes.ids
#cat insect_43_genomes.ids | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element AccessionVersion,TaxId > insect_43_genomes.taxidmapping
mmseqs createdb $DB/DB_source_fasta/insect_43_genomes.fa mmseq.genome.insectDB
mmseqs createindex mmseq.genome.insectDB tmp
mmseqs createtaxdb mmseq.genome.insectDB tmp --ncbi-tax-dump ncbi-taxdump --tax-mapping-file insect_43_genomes.taxidmapping

mmseqs createdb $DB/DB_source_fasta/insect_24_proteomes.fa mmseq.proteome.insectDB
mmseqs createindex mmseq.proteome.insectDB
