### Job name
#PBS -N LGTnoA_${id}
### Output files
#PBS -e LGTnoA_${id}.err
#PBS -o LGTnoA_${id}.log
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes/cores
#PBS -l nodes=1:ppn=1:thinnode
### Minimum memory
#PBS -l mem=10gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds>
#PBS -l walltime=4:00:00

#########################################################
# loading necessary modules                             #
#########################################################

module load tools perl samtools/1.9 bedtools/2.28.0 pigz/2.3.4 mmseqs2/release_12-113e3 barrnap/0.7 emboss/6.6.0 minimap2/2.17r941 gcc intel/perflibs R/3.6.1 lftp/4.9.2

#########################################################
# setup variables and folder structure                  #
#########################################################

# Starting time/date
STARTTIME=$(date)
STARTTIME_INSEC=$(date +%s)

# variables are passed through commandline option (-v):
# GAGA-ID: ex. -v "id=GAGA-0024"
# Squencing technology: ex. -v "tech=pacbio"

# set base directory for each genome to analyze
base=/home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${id}/

# tool directory
toolsDir=/home/projects/ku_00039/people/dinghe/github/

# set variables pointing to the final GAGA genome assembly
#assembly_dr=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Final_PacBio_assemblies_dupsrm/
#genome_file_name=$(ls -l $assembly_dr | awk -v pat="${id}" '$0~pat' | awk '{split($0,a," "); print a[9]}')
#genome=${assembly_dr}${genome_file_name}

# GAGA genome pacbio raw reads folder
#if [[ ${id} =~ ^GAGA.*$ ]]
#then
#  raw_reads_dr=/home/projects/ku_00039/people/dinghe/data/GAGA/Raw_genome_reads
#else
#  raw_reads_dr=/home/projects/ku_00039/people/dinghe/data/GAGA/Raw_genome_reads/ncbi_sra
#fi

# location of targetDB
#targetBlastnDB=/home/projects/ku_00039/people/dinghe/BLASTdb/mmseqBlastnTargetDB
#targetBlastxDB=/home/projects/ku_00039/people/dinghe/BLASTdb/swiss_prot

#################################################################################################
# switch to results folder
#################################################################################################

#base=/home/projects/ku_00039/people/luksch/GAGA/LGT/${id}/
#base=/Users/lukas/sciebo/Projects/LGT/results/${id}

cd ${base}/results
#################################################################################################

#################################################################################################
# transform blast to bed
#################################################################################################
### m6 format: query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname


#cat query.${id}.DB.euk_noAnts_blastn.m6|awk -F $'\t' 'BEGIN {OFS = FS} {if ($2<$3){print $1"\t"$2"\t"$3"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}else{print $1"\t"$3"\t"$2"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}}' > windows.fa.DB.noA.bed
tar -xOf query.${id}.DB.m6.tar.gz query.${id}.DB.euk_noAnts_blastn.m6|awk -F $'\t' 'BEGIN {OFS = FS} {if ($2<$3){print $1"\t"$2"\t"$3"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}else{print $1"\t"$3"\t"$2"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}}' > windows.fa.DB.noA.bed
#################################################################################################
# Against how many species do we get hits?
#################################################################################################

cat windows.fa.DB.noA.bed|awk -F $'\t' 'BEGIN {OFS = FS} {if ($5>69) print $1,$2,$3,$4}'|perl -pe 's/\t.*\;/\t/g'|uniq |cut -f 1|uniq -c|perl -pe 's/ +([0-9]+) (.*)/$2\t$1/g'> hitcount.species.noA.tsv

#################################################################################################
# Intersect euk and pro blast hits
#################################################################################################
# get overlapping HSPs between pro (tmp1) and euk (tmp2).
## for each pro hit, all overlapping euk hits are returned (bedtools intersect -a windows.fa.DB.pro.bed -b windows.fa.DB.euk.bed -wao)
### sort overlap by euk bitscore (sort -t$'\t' -k 1,4 -k 10,10gr)
### keep only the best (-u) scoring euk HSP (as sorted before) (-k10,10gr) for each pro HSP ( sort -t$'\t' -u -k1,4)
#################################################################################################

bedtools intersect -a windows.fa.DB.pro.bed -b windows.fa.DB.noA.bed -wao |sort -t$'\t' -k 1,4 -k 10,10gr | sort -t$'\t' -u -k1,4  > windows.fa.DB.proVSnoA.bed

#################################################################################################
# Filter 1
#################################################################################################
# Filter hsps where bacterial bitscore is higher than eukaryotic and where bacterial bitscore is > 69.
## bitscore 69 ~ evalue 1e-5
## bitscore 73 ~ evalue 1e-10
#################################################################################################

##blastn
cat windows.fa.DB.proVSnoA.bed   |awk -F $'\t' '{if ($5>$10 && $5 > 69) print $0}'|sort -t$'\t' -k 1,3 -k 5,5gr | sort -t$'\t' -u -k1,3 > windows.fa.DB.LGTs.nA.bed

#################################################################################################
# Filter 2
#################################################################################################
# Filter LGT candidate windows
## keep only those where
### scaffold is not among the "contaminants"
### the bitscore difference is > 100
### the prokaryotic HSP is larger than 100 bp
#################################################################################################

## prepare grep-able list of contaminating scaffolds
### the output in the file "contaminants" needs to be all the contaminating scaffolds in the following patterns:
#### ^scaffold1017:
#### ^scaffold1161:
#### ^scaffold1288:

### loose filter: pro HSP length > 50 bp  & bitscore diff > 50
bedtools sort -i windows.fa.DB.LGTs.nA.bed  |sort -k 1,4 -u|egrep -f contaminants - -v|awk -F $'\t' 'BEGIN {OFS = FS} {print $0"\t"$5-$10}'|sort -k 12,12gr|awk -F $'\t' 'BEGIN {OFS = FS} {if ($12 > 50 && sqrt(($2-$3)^2) > 50) print $0}' > windows.fa.DB.LGTs.nA.filtered.loose.bed

#################################################################################################
# merge HSPs to loci
#################################################################################################
# split window coordinates in tabs (perl -pe 's/(^scaffold[0-9]+)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g')
# generate genomic position (awk  -F $'\t' 'BEGIN {OFS = FS} {print $1,$2+$4,$2+$5-1,$6,$7,$12,$14,$13}')
# merge loci <500 bp apart (bedtools merge -i - -c 4,5,6,7,8 -o first,collapse,collapse,collapse,collapse -d 500)
## report first pro hit (col 4)
## collapse bac bitscores (col 5)
## collapse euk bitscores (col 6)
## collapse pairwise differences in bitscores (col 7)
## collapse overlap (col 8)
#################################################################################################

# colnames of output
# scaffold  start end bestProHit;tstart;tend;evalue;bits;alnlen;pident;taxname; collapsed(proBit)  collapsed(eukBit)  collapsed(bitDiff) collapsed(overlapping)

##blastn

cat windows.fa.DB.LGTs.nA.filtered.loose.bed| perl -pe 's/(^.+?)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g'|awk  -F $'\t' 'BEGIN {OFS = FS} {print $1,$2+$4,$2+$5-1,$6,$7,$12,$14,$13}'|bedtools sort -i - |bedtools merge -i - -c 4,5,6,7,8 -o first,collapse,collapse,collapse,collapse -d 100 |bedtools sort -i - > LGTs.nA.candidateloci.loose.bed

#################################################################################################
# Compare euk vs noAnts
#################################################################################################

bedtools intersect -a LGTs.nA.candidateloci.loose.bed -b LGTs.candidateloci.loose.bed -v > LGTs.nAo.candidateloci.loose.bed

#################################################################################################
# Analyses
#################################################################################################

# merge loci <5 kb apart
## blastn

## retrieve overlapping blastx hits
cat loci.DB.prX.bed| sort -t$'\t' -u -k1,2 |bedtools sort -i -|bedtools merge -i - -c 4 -o collapse|bedtools intersect -a LGTs.nAo.candidateloci.loose.bed -b - -wao |bedtools sort -i > LGTs.nAo.candidateloci.loose.proteins.bed

## extract regions to fasta files
#ln -s /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${id}/genome.fa ${base}/
#ln -s /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${id}/genome.fa.fai ${base}/
bedtools getfasta -fi ${base}/genome.fa -bed LGTs.nAo.candidateloci.loose.bed -fo LGTs.nAo.candidateloci.loose.fa
bedtools getfasta -fi ${base}/genome.fa -bed LGTs.nA.candidateloci.loose.bed -fo LGTs.nA.candidateloci.loose.fa

## retrieve coverage in surrounding of loci
### plus/minus 20kb, i.e. 10 windows
bedtools slop -b 20000 -i LGTs.nAo.candidateloci.loose.bed -g genome.file|bedtools intersect -b genome.overlappingwindows.cov.tsv -a - -wao |bedtools merge -i - -c 4,10,11,12,13,14 -o first,collapse,collapse,collapse,collapse,collapse|bedtools sort -i > LGTs.nAo.candidateloci.loose.coverage.bed

# Retrieve subset from bam file for each candidate
## This could be improved to produce one bam file for each candidate region.
## Should be easy with parallel, but its too late to implement now.

#BAMbase=/home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${id}/
#bedtools intersect -abam ${base}/mapping/${id}.longread.bam -b LGTs.5kb.candidateregions.bed > LGTs.nAo.5kb.candidateregions.PacBio.bam
bedtools intersect -abam ${base}/mapping/${id}.longread.bam -b LGTs.nAo.candidateloci.loose.bed > LGTs.nAo.candidateloci.loose.PacBio.bam

#################################################################################################
# ## implement screen for low complexity
#################################################################################################
# https://github.com/caballero/SeqComplex

perl -I $toolsDir/SeqComplex/ $toolsDir/SeqComplex/profileComplexSeq.pl LGTs.nAo.candidateloci.loose.fa
