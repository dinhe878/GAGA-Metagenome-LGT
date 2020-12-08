#!/bin/bash

#########################################################
# loading computerome modules
#########################################################

module load tools perl edirect/7.50 mmseqs2/release_12-113e3

#########################################################
# prepare taxidmap for mmseq database
#########################################################

# location of BLASTdb
DB=/home/people/dinghe/ku_00039/people/dinghe/BLASTdb/mmseq/

cd $DB
mkdir ncbi-taxdump && cd ncbi-taxdump
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xzvf taxdump.tar.gz
cd -

grep ">" $DB/DB_source_fasta/prok_1908_genomes.fa | awk -F ' ' '{ gsub(">","",$1); print $1 }' > prok_1908_genomes.ids
cat prok_1908_genomes.ids | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element AccessionVersion,TaxId > prok_1908_genomes.taxidmapping
mmseqs createtaxdb mmseq.genome.prokDB tmp --ncbi-tax-dump ncbi-taxdump --tax-mapping-file prok_1908_genomes.taxidmapping
