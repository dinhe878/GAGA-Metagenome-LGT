# load the compoterome modules
module load tools edirect/7.50

# go to the BLASTdb folder
cd /home/people/dinghe/ku_00039/people/dinghe/BLASTdb

# Download list of complete bacterial genome sequences
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

# generate esearch query for all chromosome level (highest assemly quality) genomes: 3307 on 20 Nov 2020
# cat prokaryotes.txt|awk -F '\t' '{if ($16=="Chromosome") print $0}'|wc -l
query=$(cat prokaryotes.txt|awk -F '\t' '{if ($16=="Chromosome") print $0}'| |cut -f 19|sed 1d|tr "\n" " "|perl -pe 's/ / OR /g'|perl -pe "s/ OR $//g")

# retrieve genomes
esearch -db assembly -query "$query"|elink -target nuccore|efetch -format fasta > db.prok_3307genomes_20112020.fa

# retrieve proteomes
esearch -db assembly -query "$query"|elink -target nuccore|efetch -format fasta_cds_aa > db.prok_3307genomes_20112020.proteins.fa
