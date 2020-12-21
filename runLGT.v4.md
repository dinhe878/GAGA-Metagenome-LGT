# Detecting LGTs and prokaryotic scaffolds in draft genomes

<!-- TOC START min:1 max:3 link:true asterisk:false update:true -->
- [Detecting LGTs and prokaryotic scaffolds in draft genomes](#detecting-lgts-and-prokaryotic-scaffolds-in-draft-genomes)
- [Dependencies](#dependencies)
  - [Setup variables](#setup-variables)
  - [Prepare folder structure](#prepare-folder-structure)
- [Setup bin folder](#setup-bin-folder)
  - [Prepare input files](#prepare-input-files)
  - [Split genome in windows to process](#split-genome-in-windows-to-process)
    - [Run samtools faidx to retrieve windows](#run-samtools-faidx-to-retrieve-windows)
  - [mmseq blastn vs a prokaryote DB and an insect DB](#mmseq-blastn-vs-a-prokaryote-db-and-an-insect-db)
    - [Create mmseqs dbs](#create-mmseqs-dbs)
    - [Run mmseq](#run-mmseq)
- [rRNA prediction](#rrna-prediction)
  - [Calculate GC content and length for each scaffold](#calculate-gc-content-and-length-for-each-scaffold)
- [Run mmseq blastx vs protein dbs](#run-mmseq-blastx-vs-protein-dbs)
- [Analyze long-read coverage](#analyze-long-read-coverage)
- [Retrieve results files](#retrieve-results-files)
- [Analyze results](#analyze-results)
- [Supplement](#supplement)
  - [Downloading sets of prokaryotic genomes and proteomes (n=2422)](#downloading-sets-of-prokaryotic-genomes-and-proteomes-n2422)
    - [Data download](#data-download)
    - [Database preparation](#database-preparation)
  - [Downloading sets of HQ insect genomes and proteomes](#downloading-sets-of-hq-insect-genomes-and-proteomes)
    - [Data download](#data-download-1)
    - [Database preparation](#database-preparation-1)
<!-- TOC END -->

# Dependencies
```bash
samtools v1.10
bedtools
barrnap
mmseqs
infoseq
minimap2

eutils #for retrieving genomes and proteomes

```


**This workflow will:**
- blast genomic overlapping windows of a certain size (e.g. 3.5 kb) against a
  - **prokaryotic genome database with mmseqs blastn**
  - **eukaryotic genome database with mmseqs blastn**
  - **prokaryotic protein database with mmseqs blastx**
  - **insect protein database with DIAMOND blastX**

- annotated **eukaryotic rRNAs** and **prokaryotic rRNAs** with barrnap



## Setup variables
```bash
id=GAGA-0361
# set base directory for each genome to analyze
base=/home/people/dinghe/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/$id/

# GAGA genome is the polished assembly
genome=$(readlink -f /home/people/dinghe/ku_00039/people/dinghe/data/GAGA/Polished_assemblies/$id*_nextpolish.fasta.gz)

# location of BLASTdb
DB=/home/people/dinghe/ku_00039/people/dinghe/BLASTdb/mmseq/
```

## Prepare folder structure
```bash
mkdir $base
cd $base
mkdir $base/mapping/
mkdir $base/results
mkdir $base/mmseqs
```

## Prepare input files
```bash
# prepare genome fasta
zcat $genome|cut -f 1 -d "|"  > $base/genome.fa
samtools faidx genome.fa
```
## Split genome in windows to process

```bash
# generate 2.5 kb windows with 500 bp walking steps (overlap 500 bp with the previous and the subsequent window)

# calculate windows
cut -f 1,2 genome.fa.fai > genome.lengths.tsv
windowSize=2000
stepSize=1500

# generate overlapping windows region file (for samtools faidx)
./bin/bedtools makewindows -g genome.lengths.tsv -w $windowSize -s $stepSize > genome.overlappingwindows.bed
cat genome.overlappingwindows.bed|perl -pe 's/(.*?)\t(.*?)\t(.*?)/$1\:$2\-$3/g'  > genome.overlappingwindows.tsv


```

### Run samtools faidx to retrieve windows
Here, we create a fasta file for all the windows we will process.
```bash
# Create a single fasta file with all windows
#requries samtools v1.10
./bin/samtools faidx genome.fa --region genome.overlappingwindows.tsv > windows.fa
```

The file `windows.fa` contains all overlapping windows.

## mmseq blastn vs a prokaryote DB and an insect DB
### Create mmseqs dbs
```bash
cd $base

# index windows.fa for use with mmseqs
./bin/mmseqs createdb windows.fa windows.fa.DB
```
### Run mmseq
Now, **all overlapping windows** will be blasted against the prokaryotic database and then the eukaryotic databse with mmseqs. Hits against the prokaryotic db are possible LGTs and bacterial scaffolds/contigs.

This is done by running `mmseqs search` locally.

```bash
# define windows fasta mmseqs database
inDB=windows.fa.DB

# define mmseqs search sensitivity
sensitivity=7

# set $tag and $db for the prokaryotic screen
tag="pro"
db="/global/scratch/schradel/BLAST/dbs/mmseq.prokDB"

# set $tag and $db for the eukaryotic screen
tag="euk"
db="/global/scratch/schradel/BLAST/dbs/mmseq.insectDB"

# run mmseqs search (twice, first on prokaryotic (i.e. tag="pro") and then on eukaryotic (tag="euk"))
# control resource usage with --threads and --split-memory-limit when running on a cluster
time ./bin/mmseqs search ${base}/${inDB} ${db} $base/mmseqs/${inDB}.${tag}.resDB tmp --start-sens 1 --sens-steps 2 -s ${sensitivity} --search-type 3 1> $base/mmseqs/${inDB}.${tag}.search.out 2> $base/mmseqs/${inDB}.${tag}.search.err

# run mmseqs convertalis (to generate blast m8-like  output format)
time ./bin/mmseqs convertalis ${base}/${inDB} ${db} $base/mmseqs/${inDB}.${tag}.resDB $base/mmseqs/${inDB}.${tag}.m6 --format-output query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname  1> $base/mmseqs/${inDB}.${tag}.convert.out 2> $base/mmseqs/${inDB}.${tag}.convert.err

# sort first by evalue (-k7,7n), then by bitscore (-k8,8nr)
# keep best hit only
cat $base/mmseqs/${inDB}.${tag}.m6 |sort -k1,1 -k7,7n -k8,8nr | sort -u -k1,1 --merge |perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2n > $base/mmseqs/${inDB}.${tag}.bh

# retrieve a list of all the windows that hit against the prokDB
cat $base/mmseqs/${inDB}.${tag}.bh |cut -f 4 > $base/mmseqs/${inDB}.${tag}.bh.lst

```

# rRNA prediction
This is done with **barrnap** using the [batch submission script `batchBlast.rRNA.qsub.sh`](#batch-script-for-rrna-prediction-both-eukaryots-and-prokaryots) .
Run this batch script to annotate rRNA in those windows that hit against the prokDB. This is necessary because many of the hits against the prokaryotic DB are in fact eukaryotic rRNAs. These are so universally conserved that they also hit against the prokaryotic genomes.

```bash
cd $base
# prokaryotic and eukaryotic rRNAs
## annotate rRNAs
barrnap --threads 20 -k bac genome.fa > genome.pro.rRNA.gff3
barrnap --threads 20 -k euk genome.fa > genome.euk.rRNA.gff3
## retrieve windows that overlap an rRNA
bedtools intersect -a genome.windows.tsv -b genome.pro.rRNA.gff3 -wa -wb > genome.pro.rRNA.windows.bed
bedtools intersect -a genome.windows.tsv -b genome.euk.rRNA.gff3 -wa -wb > genome.euk.rRNA.windows.bed
```

## Calculate GC content and length for each scaffold
Calculate **GC content** for all scaffolds and for all windows for the subsequent analysis in R.
```bash
# calculate gc and length by scaffold
infoseq  -nocolumn -delimiter "\t" -auto -only  -name -length -pgc genome.fa > $base/genome.GC.tsv

# calculate gc and length by window
cat windows.fa|perl -pe 's/(>.*?)\:(.*)/$1\@$2/g' > tmp.fa
infoseq  -nocolumn -delimiter "\t" -auto -only  -name -length -pgc tmp.fa |perl -pe 's/^(.*?)@(.*)/$1:$2/g' > ./windows.GC.tsv
rm tmp.fa
```

# Run mmseq blastx vs protein dbs
All windows will be compared against a prokaryot protein database with mmseqs blastx.
```bash
cd $base
#define windows mmseqs database
inDB=windows.fa.DB

# define mmseqs search sensitivity
sensitivity=7

# set $tag and $db for the prokaryotic PROTEIN screen
tag="proP"
db="/global/scratch/schradel/BLAST/dbs/mmseq.prokPDB"

# set $tag and $db for the eukaryotic PROTEIN screen
tag="eukP"
db="/global/scratch/schradel/BLAST/dbs/mmseq.insectPDB"

# run mmseqs search
# control resource usage with --threads and --split-memory-limit when running on a cluster
time ${mmseqBin} search ${base}/${inDB} ${db} $base/mmseqs/${inDB}.${tag}.resDB tmp --start-sens 1 --sens-steps 2 -s ${sensitivity} --search-type 3 1> $base/mmseqs/${inDB}.${tag}.search.out 2> $base/mmseqs/${inDB}.${tag}.search.err

# run mmseqs convertalis (to generate blast m8-like  output format)
time ${mmseqBin} convertalis ${base}/${inDB} ${db} $base/mmseqs/${inDB}.${tag}.resDB $base/mmseqs/${inDB}.${tag}.m6 --format-output query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident  1> $base/mmseqs/${inDB}.${tag}.convert.out 2> $base/mmseqs/${inDB}.${tag}.convert.err

# select best blastx hit per window
cat $base/mmseqs/${inDB}.${tag}.m6|sort -k1,1 -k7,7n -k8,8nr | sort -u -k1,1 --merge |perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2n > $base/mmseqs/${inDB}.${tag}.bh

```


# Analyze long-read coverage
```bash

cd $base/mapping

# retrieve long reads
ln -s $longRead.fq.gz ./${id}.raw.fq.gz

cd $base
# map to genome with minimap2
./bin/minimap2 -t 5 -ax map-ont $base/genome.fa $base/mapping/${id}.raw.fq.gz > $base/mapping/${id}.longread.sam

#sort sam and convert to bam
./bin/samtools view -S -b $base/mapping/${id}.longread.sam |./bin/samtools sort > $base/mapping/${id}.longread.bam

# calculate coverage for each window
coverageBed -a $base/genome.overlappingwindows.bed -b $base/mapping/${id}.longread.bam > $base/mapping/genome.overlappingwindows.cov.tsv
```

# Retrieve results files
All relevant results are finally stored in ```$base/results/```
```
ln -s $base/genome.GC.tsv results/
ln -s $base/windows.GC.tsv results/
# rRNA predictions
ln -s $base/genome.pro.rRNA.windows.bed results/
ln -s $base/genome.euk.rRNA.windows.bed results/
# window blastn
ln -s $base/mmseqs/windows.fa.DB.pro.bh results/
ln -s $base/mmseqs/windows.fa.DB.euk.bh results/
ln -s $base/genome.windows.tsv results/
ln -s $base/mapping/genome.overlappingwindows.cov.tsv results/

# protein blastx
ln -s $base/mmseqs/windows.fa.DB.proP.bh results/
ln -s $base/mmseqs/windows.fa.DB.eukP.bh results/

```

# Analyze results
These results can be analyzed in R to identify contaminating scaffolds and to identify putative LGTs.
**[see ```analyseBAC.Rmd```](./analyseBAC.Rmd)**


# Supplement
## Downloading sets of prokaryotic genomes and proteomes (n=2422)
### Data download
Download a set of prokaryotic genomes+proteomes that are assembled to Chromosome level. Also retrieve only one genome per subgroup.

```bash
# download proteins for all 2422 prok. genomes
cd /global/scratch/schradel/BLAST/dbs
# Download list of complete bacterial genome sequences
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

# retrieve unique subgroup (field 7)(or group (field 6)) for taxonomic diversity
cat /global/scratch/schradel/BLAST/dbs/prokaryotes.txt|awk -F '\t' '{if ($16=="Chromosome") print $0}'|awk -F '\t' '!seen[$7]++'|wc -l

# % 2422 genomes
query=$(cat /global/scratch/schradel/BLAST/dbs/prokaryotes.txt|awk -F '\t' '{if ($16=="Chromosome") print $0}'|awk -F '\t' '!seen[$7]++'|cut -f 19|sed 1d|tr "\n" " "|perl -pe 's/ / OR /g'|perl -pe "s/ OR $//g")

# retrieve genomes
esearch -db assembly -query "$query"|elink -target nuccore|efetch -format fasta > db.2422genomes.fa

# retrieve proteomes
esearch -db assembly -query "$query"|elink -target nuccore|efetch -format fasta_cds_aa > db.2422genomes.proteins.fa
```

### Database preparation
```bash
# format db
cd /global/scratch/schradel/BLAST/dbs/

# create mmseqs db from prokaryotic genomes fasta file
mmseqs createdb db.2422genomes.fa mmseq.prokDB

## don't create an mmseqs index or the SGE runs will fail when run on a cluster
# mmseqs createindex mmseq.prokDB tmp --search-type 3 --split-memory-limit 30G

# create diamond blastx DB for proteins
mmseqs createdb db.2422genomes.proteins.fa mmseq.prokPDB

# prepare taxidmap for mmseq database
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
mkdir new_taxonomy && tar -xxvf new_taxdump.tar.gz -C new_taxonomy

# extract ids
grep ">" db.2422genomes.fa |perl -pe 's/\>(.*?) .*/$1/g' > db.2422genomes.ids
# retrieve taxonomy with eutils
cat db.2422genomes.ids| epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element AccessionVersion,TaxId > db.2422genomes.taxidmap

# create taxDB for mmseqs db
rm mmseq.prokDB_mapping
mmseqs createtaxdb mmseq.prokDB tmp --ncbi-tax-dump new_taxonomy --tax-mapping-file db.2422genomes.taxidmap

# the file should be larger than 0 bytes
ls -lh mmseq.prokDB_mapping


```
## Downloading sets of HQ insect genomes and proteomes
### Data download
Download genome+proteomes for a set of high quality insect genome assemblies with status either "Chromosome" or "Complete genome".
```bash
# get list of all genomes from NCBI (~30.11.2019)
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt
# get unique insect genera that were fully assembled
query=$(cat /global/scratch/schradel/BLAST/dbs/eukaryotes.txt|awk -F '\t' '{if ($17=="Chromosome" || $17=="Complete Genome") print $0}'|awk -F ' ' '!seen[$1]++'|awk -F '\t' '{if ($6=="Insects") print $0}'|cut -f 9|sed 1d|tr "\n" " "|perl -pe 's/ / OR /g'|perl -pe "s/ OR $//g")

# retrieve genomes
esearch -db assembly -query "$query"|elink -target nuccore|efetch -format fasta > db.insectGenomes.fa

# retrieve fasta of proteomes
esearch -db assembly  -query "$query"|elink -target nuccore|efetch -format fasta_cds_aa > db.insectGenomes.proteins.fa
```

### Database preparation
```bash
# make blast DB
cd /global/scratch/schradel/BLAST/dbs/

mmseqs createdb db.insectGenomes.fa mmseq.insectDB

## don't create an mmseqs index or the SGE runs will fail when run on a cluster
# mmseqs createindex mmseq.insectDB tmp2 --search-type 3 --split-memory-limit 30G

# create diamond blastx DB for proteins
mmseqs createdb db.insectGenomes.proteins.fa mmseq.insectPDB

# prepare taxidmap for mmseq database
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
mkdir new_taxonomy && tar -xxvf new_taxdump.tar.gz -C new_taxonomy

# extract ids
grep ">" db.insectGenomes.fa |perl -pe 's/\>(.*?) .*/$1/g' > db.insectGenomes.ids
# retrieve taxonomy with eutils (split into multiple jobs with parallel)
export NCBI_API_KEY='7820ce7ebdfc0ed36d2d1d2fc47746f0c209'
cat db.insectGenomes.ids| parallel --delay 10 --pipe -N 1000 -j 3 "cat|epost -db nuccore |esummary -db nuccore |xtract -pattern DocumentSummary -element AccessionVersion,TaxId" > db.insectGenomes.taxidmap

# create taxDB for mmseqs db
rm mmseq.insectDB_mapping
mmseqs createtaxdb mmseq.insectDB tmp --ncbi-tax-dump new_taxonomy --tax-mapping-file db.insectGenomes.taxidmap

# the file should be larger than 0 bytes
ls -lh mmseq.insectDB_mapping

```
