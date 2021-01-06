### Job name
#PBS -N LGT_GAGA-0090
### Output files
#PBS -e LGT_GAGA-0090.err
#PBS -o LGT_GAGA-0090.log
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes/cores
#PBS -l nodes=1:ppn=40:thinnode
### Minimum memory
#PBS -l mem=150gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds>
#PBS -l walltime=12:00:00

#########################################################
# loading necessary modules                             #
#########################################################

module load tools perl samtools/1.10 bedtools/2.28.0 pigz/2.3.4 mmseqs2/release_12-113e3 barrnap/0.7 emboss/6.6.0 minimap2/2.17r941

#########################################################
# setup variables and folder structure                  #
#########################################################

# Starting time/date
STARTTIME=$(date)
STARTTIME_INSEC=$(date +%s)

# GAGA-ID
id=GAGA-0090

# set base directory for each genome to analyze
base=/home/people/dinghe/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${id}/

# GAGA genome is the polished assembly
genome=$(readlink -f /home/people/dinghe/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Verified_polished_assemblies/${id}_nextpolish.fasta)

# GAGA genome pacbio raw reads folder
raw_reads_dr=/home/people/dinghe/ku_00039/people/dinghe/data/GAGA/Raw_genome_reads

# location of batch scripts
scripts=/home/people/dinghe/ku_00039/people/dinghe/scripts/batch/

# location of BLASTdb
DB=/home/people/dinghe/ku_00039/people/dinghe/BLASTdb/mmseq/

mkdir $base
cd $base
mkdir mapping
mkdir results
mkdir mmseqs
mkdir tmp        # if not already exist

#########################################################
# prepare input files                                   #
#########################################################

# prepare genome fasta
cat $genome|cut -f 1 -d "|"  > genome.fa
samtools faidx genome.fa

# generate 2.5 kb windows with 500 bp walking steps (overlap 500 bp with the previous and the subsequent window)
# calculate windows
cut -f 1,2 genome.fa.fai > genome.lengths.tsv
windowSize=2000
stepSize=1500
bedtools makewindows -g genome.lengths.tsv -w $windowSize -s $stepSize > genome.overlappingwindows.bed
cat genome.overlappingwindows.bed|perl -pe 's/(.*?)\t(.*?)\t(.*?)/$1\:$2\-$3/g'  > genome.overlappingwindows.tsv
samtools faidx genome.fa --region genome.overlappingwindows.tsv > windows.fa

# create mmseqs query database
mmseqs createdb windows.fa query.${id}.DB

# define windows fasta mmseqs database (query)
inDB=query.${id}.DB

# define mmseqs search sensitivity
sensitivity=7

# set tag and db for the prokaryotic screen
tag_pro="pro"
db_pro="/home/people/dinghe/ku_00039/people/dinghe/BLASTdb/mmseq/mmseq.genome.prokDB"

# set tag and db for the eukaryotic screen
tag_euk="euk"
db_euk="/home/people/dinghe/ku_00039/people/dinghe/BLASTdb/mmseq/mmseq.genome.insectDB"

# run mmseqs search twice, first on prokaryotic (i.e. tag="pro") and then on eukaryotic (tag="euk")
# if applicable, control resource usage with --threads and --split-memory-limit when running on HPC
mmseqs search ${inDB} ${db_pro} mmseqs/${inDB}.${tag_pro}.resDB tmp --start-sens 1 --sens-steps 2 -s ${sensitivity} --search-type 3 1> mmseqs/${inDB}.${tag_pro}.search.out 2> mmseqs/${inDB}.${tag_pro}.search.err
mmseqs search ${inDB} ${db_euk} mmseqs/${inDB}.${tag_euk}.resDB tmp --start-sens 1 --sens-steps 2 -s ${sensitivity} --search-type 3 1> mmseqs/${inDB}.${tag_euk}.search.out 2> mmseqs/${inDB}.${tag_euk}.search.err

# run mmseqs convertalis (to generate blast m6-like output format)
mmseqs convertalis ${inDB} ${db_pro} mmseqs/${inDB}.${tag_pro}.resDB mmseqs/${inDB}.${tag_pro}.m6 --format-output query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname  1> mmseqs/${inDB}.${tag_pro}.convert.out 2> mmseqs/${inDB}.${tag_pro}.convert.err

mmseqs convertalis ${inDB} ${db_euk} mmseqs/${inDB}.${tag_euk}.resDB mmseqs/${inDB}.${tag_euk}.m6 --format-output query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname  1> mmseqs/${inDB}.${tag_euk}.convert.out 2> mmseqs/${inDB}.${tag_euk}.convert.err

# sort first by evalue (-k7,7g), then by bitscore (-k8,8gr)
# keep best hit only
cat mmseqs/${inDB}.${tag_pro}.m6 |sort -k1,1 -k7,7g -k8,8gr | sort -u -k1,1 --merge |perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > mmseqs/${inDB}.${tag_pro}.bh
cat mmseqs/${inDB}.${tag_euk}.m6 |sort -k1,1 -k7,7g -k8,8gr | sort -u -k1,1 --merge |perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > mmseqs/${inDB}.${tag_euk}.bh

# retrieve a list of all the windows that hit against the prokDB
cat mmseqs/${inDB}.${tag_pro}.bh |cut -f 4 > mmseqs/${inDB}.${tag_pro}.bh.lst
cat mmseqs/${inDB}.${tag_euk}.bh |cut -f 4 > mmseqs/${inDB}.${tag_euk}.bh.lst

# rRNA prediction
barrnap --threads 40 -k bac genome.fa > genome.${tag_pro}.rRNA.gff3
barrnap --threads 40 -k euk genome.fa > genome.${tag_euk}.rRNA.gff3
## retrieve windows that overlap an rRNA
bedtools intersect -a genome.overlappingwindows.bed -b genome.${tag_pro}.rRNA.gff3 -wa -wb > genome.${tag_pro}.rRNA.windows.bed
bedtools intersect -a genome.overlappingwindows.bed -b genome.${tag_euk}.rRNA.gff3 -wa -wb > genome.${tag_euk}.rRNA.windows.bed

# Calculate GC content and length for each scaffold
infoseq  -nocolumn -delimiter "\t" -auto -only -name -length -pgc genome.fa > genome.GC.tsv

# calculate gc and length by window
cat windows.fa|perl -pe 's/(>.*?)\:(.*)/$1\@$2/g' > tmp.fa
infoseq  -nocolumn -delimiter "\t" -auto -only -name -length -pgc tmp.fa |perl -pe 's/^(.*?)@(.*)/$1:$2/g' > ./windows.GC.tsv
rm tmp.fa

# Analyze long-read coverage
bamToFastq -i ${raw_reads_dr}/bam/${id}.bam -fq ${raw_reads_dr}/fq/${id}.fq.gz
minimap2 -t 40 -ax map-pb genome.fa ${raw_reads_dr}/fq/${id}.fq.gz > mapping/${id}.longread.sam
samtools view -S -b mapping/${id}.longread.sam | samtools sort > mapping/${id}.longread.bam
bedtools coverage -a genome.overlappingwindows.bed -b mapping/${id}.longread.bam > mapping/genome.overlappingwindows.cov.tsv

# Gather results
# GC content
mv genome.GC.tsv results/
mv windows.GC.tsv results/
# rRNA predictions
mv genome.pro.rRNA.windows.bed results/
mv genome.euk.rRNA.windows.bed results/
# window blastn
mv mmseqs/${inDB}.${tag_pro}.bh results/
mv mmseqs/${inDB}.${tag_euk}.bh results/
mv genome.overlappingwindows.bed results/
mv mapping/genome.overlappingwindows.cov.tsv results/

# clean up tmp files to save space
rm -rf tmp/*

# Ending time/date
ENDTIME=$(date)
ENDTIME_INSEC=$(date +%s)
echo "==============================================="
echo "Pipeline started at $STARTTIME"
echo "Pipeline ended at $ENDTIME"
echo "Pipeline took $((ENDTIME_INSEC - STARTTIME_INSEC)) seconds to finish"
