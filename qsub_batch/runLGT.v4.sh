### Job name
#PBS -N Bac_screen_${id}
### Output files
#PBS -e Bac_screen_${id}.err
#PBS -o Bac_screen_${id}.log
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes/cores
#PBS -l nodes=1:ppn=40:thinnode
### Minimum memory
#PBS -l mem=180gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds>
#PBS -l walltime=24:00:00

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
# Check if only assembled with stLFR: -v "stLFR=T"

# set base directory for each genome to analyze
base=/home/people/dinghe/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${id}/

# tool directory
toolsDir=/home/people/dinghe/ku_00039/people/dinghe/github/

# set variables pointing to the final GAGA genome assembly
if [[ ${stLFR} == T ]]
then
  assembly_dr=/home/people/dinghe/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Final_stLFR_assemblies_dupsrm/
else
  assembly_dr=/home/people/dinghe/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Final_PacBio_assemblies_dupsrm/
fi

genome_file_name=$(ls -l $assembly_dr | awk -v pat="${id}" '$0~pat' | awk '{split($0,a," "); print a[9]}')
genome=${assembly_dr}${genome_file_name}

# GAGA genome pacbio raw reads folder
if [[ ${id} =~ ^GAGA.*$ ]]
then
  raw_reads_dr=/home/people/dinghe/ku_00039/people/dinghe/data/GAGA/Raw_genome_reads
else
  raw_reads_dr=/home/people/dinghe/ku_00039/people/dinghe/data/GAGA/Raw_genome_reads/ncbi_sra
fi

# location of targetDB
targetBlastnDB=/home/people/dinghe/ku_00039/people/dinghe/BLASTdb/mmseqBlastnTargetDB
targetBlastxDB=/home/people/dinghe/ku_00039/people/dinghe/BLASTdb/swiss_prot
mkdir $base
cd $base
mkdir mapping
mkdir results
mkdir mmseqs
mkdir tmp        # if not already exist

#########################################################
# prepare input files                                   #
#########################################################

# prepare genome query fasta
echo "Preparing genome query fasta..."
cat $genome | cut -f 1 -d "|" | sed 's/:/_/g' > genome.fa
samtools faidx genome.fa

# generate 2.5 kb windows with 500 bp walking steps (overlap 500 bp with the previous and the subsequent window)
# calculate windows
cut -f 1,2 genome.fa.fai > genome.lengths.tsv
windowSize=2000
stepSize=1500
bedtools makewindows -g genome.lengths.tsv -w $windowSize -s $stepSize > genome.overlappingwindows.bed
bedtools sort -i genome.overlappingwindows.bed > genome.overlappingwindows.sorted.bed
cat genome.overlappingwindows.sorted.bed|perl -pe 's/(.*?)\t(.*?)\t(.*?)/$1\:$2\-$3/g'  > genome.overlappingwindows.tsv
samtools faidx genome.fa --region genome.overlappingwindows.tsv > windows.fa

# create mmseqs query database
echo "Creating mmseqs query database..."
mmseqs createdb windows.fa query.${id}.DB

# define windows fasta mmseqs database (query)
inDB=query.${id}.DB

# define mmseqs search sensitivity
sensitivity=7

# set tag and db for the prokaryotic screen
tag_pro_n="pro_blastn"
tag_pro_x="pro_blastx"
genome_db_pro="$targetBlastnDB/mmseq.genome.prokDB"
proteome_db_pro="$targetBlastxDB/Bacteria"

# set tag and db for the eukaryotic screen
tag_euk_n="euk_blastn"
tag_euk_x="euk_blastx"
tag_euk_noAnts_n="euk_noAnts_blastn"
tag_human_n="human_blastn"
genome_db_insect="$targetBlastnDB/mmseq.genome.71_clean_insectDB"
genome_db_noAnts_insect="$targetBlastnDB/mmseq.genome.clean.noAnts_insectDB"
genome_db_human="$targetBlastnDB/mmseq.genome.humanDB"
proteome_db_euk="$targetBlastxDB/Insecta"

# mmseqs search on bacterial db, selected insects db, selected insects db excluding ants, and human db
echo "Starting mmseqs blasting ${inDB} against ${genome_db_pro}..."
mmseqs search ${inDB} ${genome_db_pro} mmseqs/${inDB}.${tag_pro_n}.resDB tmp --start-sens 1 --sens-steps 2 -s ${sensitivity} --search-type 3 1> mmseqs/${inDB}.${tag_pro_n}.search.out 2> mmseqs/${inDB}.${tag_pro_n}.search.err
echo "Starting mmseqs blasting ${inDB} against ${proteome_db_pro}..."
mmseqs search ${inDB} ${proteome_db_pro} mmseqs/${inDB}.${tag_pro_x}.resDB tmp --start-sens 1 --sens-steps 2 -s ${sensitivity} --search-type 3 1> mmseqs/${inDB}.${tag_pro_x}.search.out 2> mmseqs/${inDB}.${tag_pro_x}.search.err
echo "Starting mmseqs blasting ${inDB} against ${genome_db_insect}..."
mmseqs search ${inDB} ${genome_db_insect} mmseqs/${inDB}.${tag_euk_n}.resDB tmp --start-sens 1 --sens-steps 2 -s ${sensitivity} --search-type 3 1> mmseqs/${inDB}.${tag_euk_n}.search.out 2> mmseqs/${inDB}.${tag_euk_n}.search.err
echo "Starting mmseqs blasting ${inDB} against ${proteome_db_euk}..."
mmseqs search ${inDB} ${proteome_db_euk} mmseqs/${inDB}.${tag_euk_x}.resDB tmp --start-sens 1 --sens-steps 2 -s ${sensitivity} --search-type 3 1> mmseqs/${inDB}.${tag_euk_x}.search.out 2> mmseqs/${inDB}.${tag_euk_x}.search.err
echo "Starting mmseqs blasting ${inDB} against ${genome_db_noAnts_insect}..."
mmseqs search ${inDB} ${genome_db_noAnts_insect} mmseqs/${inDB}.${tag_euk_noAnts_n}.resDB tmp --start-sens 1 --sens-steps 2 -s ${sensitivity} --search-type 3 1> mmseqs/${inDB}.${tag_euk_noAnts_n}.search.out 2> mmseqs/${inDB}.${tag_euk_noAnts_n}.search.err
echo "Starting mmseqs blasting ${inDB} against ${genome_db_human}..."
mmseqs search ${inDB} ${genome_db_human} mmseqs/${inDB}.${tag_human_n}.resDB tmp --start-sens 1 --sens-steps 2 -s ${sensitivity} --search-type 3 1> mmseqs/${inDB}.${tag_human_n}.search.out 2> mmseqs/${inDB}.${tag_human_n}.search.err

# run mmseqs convertalis (to generate blast m6-like output format)
echo "Converting blast results..."
mmseqs convertalis ${inDB} ${genome_db_pro} mmseqs/${inDB}.${tag_pro_n}.resDB mmseqs/${inDB}.${tag_pro_n}.m6 --format-output query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname  1> mmseqs/${inDB}.${tag_pro_n}.convert.out 2> mmseqs/${inDB}.${tag_pro_n}.convert.err
mmseqs convertalis ${inDB} ${proteome_db_pro} mmseqs/${inDB}.${tag_pro_x}.resDB mmseqs/${inDB}.${tag_pro_x}.m6 --format-output query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname  1> mmseqs/${inDB}.${tag_pro_x}.convert.out 2> mmseqs/${inDB}.${tag_pro_x}.convert.err

mmseqs convertalis ${inDB} ${genome_db_insect} mmseqs/${inDB}.${tag_euk_n}.resDB mmseqs/${inDB}.${tag_euk_n}.m6 --format-output query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname  1> mmseqs/${inDB}.${tag_euk_n}.convert.out 2> mmseqs/${inDB}.${tag_euk_n}.convert.err
mmseqs convertalis ${inDB} ${proteome_db_euk} mmseqs/${inDB}.${tag_euk_x}.resDB mmseqs/${inDB}.${tag_euk_x}.m6 --format-output query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname  1> mmseqs/${inDB}.${tag_euk_x}.convert.out 2> mmseqs/${inDB}.${tag_euk_x}.convert.err

mmseqs convertalis ${inDB} ${genome_db_noAnts_insect} mmseqs/${inDB}.${tag_euk_noAnts_n}.resDB mmseqs/${inDB}.${tag_euk_noAnts_n}.m6 --format-output query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname  1> mmseqs/${inDB}.${tag_euk_noAnts_n}.convert.out 2> mmseqs/${inDB}.${tag_euk_noAnts_n}.convert.err

mmseqs convertalis ${inDB} ${genome_db_human} mmseqs/${inDB}.${tag_human_n}.resDB mmseqs/${inDB}.${tag_human_n}.m6 --format-output query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname  1> mmseqs/${inDB}.${tag_human_n}.convert.out 2> mmseqs/${inDB}.${tag_human_n}.convert.err

# sort first by evalue (-k7,7g), then by bitscore (-k8,8gr)
# keep best hit only
cat mmseqs/${inDB}.${tag_pro_n}.m6 |sort -k1,1 -k7,7g -k8,8gr | sort -u -k1,1 --merge |perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > mmseqs/${inDB}.${tag_pro_n}.bh
cat mmseqs/${inDB}.${tag_pro_x}.m6 |sort -k1,1 -k7,7g -k8,8gr | sort -u -k1,1 --merge |perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > mmseqs/${inDB}.${tag_pro_x}.bh
cat mmseqs/${inDB}.${tag_euk_n}.m6 |sort -k1,1 -k7,7g -k8,8gr | sort -u -k1,1 --merge |perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > mmseqs/${inDB}.${tag_euk_n}.bh
cat mmseqs/${inDB}.${tag_euk_x}.m6 |sort -k1,1 -k7,7g -k8,8gr | sort -u -k1,1 --merge |perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > mmseqs/${inDB}.${tag_euk_x}.bh
cat mmseqs/${inDB}.${tag_euk_noAnts_n}.m6 |sort -k1,1 -k7,7g -k8,8gr | sort -u -k1,1 --merge |perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > mmseqs/${inDB}.${tag_euk_noAnts_n}.bh
cat mmseqs/${inDB}.${tag_human_n}.m6 |sort -k1,1 -k7,7g -k8,8gr | sort -u -k1,1 --merge |perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > mmseqs/${inDB}.${tag_human_n}.bh

# retrieve a list of all the windows that hit against the prokDB
cat mmseqs/${inDB}.${tag_pro_n}.bh |cut -f 4 > mmseqs/${inDB}.${tag_pro_n}.bh.lst
cat mmseqs/${inDB}.${tag_pro_x}.bh |cut -f 4 > mmseqs/${inDB}.${tag_pro_x}.bh.lst
cat mmseqs/${inDB}.${tag_euk_n}.bh |cut -f 4 > mmseqs/${inDB}.${tag_euk_n}.bh.lst
cat mmseqs/${inDB}.${tag_euk_x}.bh |cut -f 4 > mmseqs/${inDB}.${tag_euk_x}.bh.lst
cat mmseqs/${inDB}.${tag_euk_noAnts_n}.bh |cut -f 4 > mmseqs/${inDB}.${tag_euk_noAnts_n}.bh.lst
cat mmseqs/${inDB}.${tag_human_n}.bh |cut -f 4 > mmseqs/${inDB}.${tag_human_n}.bh.lst

# rRNA prediction
echo "Starting rRNA prediction..."
barrnap --threads 40 -k bac genome.fa > genome.${tag_pro_n}.rRNA.gff3
barrnap --threads 40 -k euk genome.fa > genome.${tag_euk_n}.rRNA.gff3
## retrieve windows that overlap an rRNA
bedtools intersect -a genome.overlappingwindows.sorted.bed -b genome.${tag_pro_n}.rRNA.gff3 -wa -wb > genome.${tag_pro_n}.rRNA.windows.bed
bedtools intersect -a genome.overlappingwindows.sorted.bed -b genome.${tag_euk_n}.rRNA.gff3 -wa -wb > genome.${tag_euk_n}.rRNA.windows.bed

# Calculate GC content and length for each scaffold
echo "Calculating GC content of scaffolds..."
infoseq  -nocolumn -delimiter "\t" -auto -only -name -length -pgc genome.fa > genome.GC.tsv

# calculate gc and length by window
cat windows.fa|perl -pe 's/(>.*?)\:(.*)/$1\@$2/g' > tmp.fa
infoseq  -nocolumn -delimiter "\t" -auto -only -name -length -pgc tmp.fa |perl -pe 's/^(.*?)@(.*)/$1:$2/g' > ./windows.GC.tsv
rm tmp.fa

# Analyze long-read coverage
echo "Gathering sequencing coverage information..."
if [[ ${id} =~ ^GAGA.*$ ]]
then
  minimap2 -t 40 -ax map-pb genome.fa ${raw_reads_dr}/fq/${id}.fq.gz > mapping/${id}.longread.sam
elif [[ ${tech} == "pair" ]]
then
  minimap2 -t 40 -ax sr genome.fa ${raw_reads_dr}/fq/${id}_1.fq.gz ${raw_reads_dr}/fq/${id}_2.fq.gz > mapping/${id}.longread.sam
elif [[ ((${id} =~ ^NCBI.*$) || (${id} =~ ^OUT.*$)) && (${tech} == "pacbio") ]]
then
  minimap2 -t 40 -ax map-pb genome.fa ${raw_reads_dr}/fq/${id}.pacbio.fq.gz > mapping/${id}.longread.sam
else
  minimap2 -t 40 -ax sr genome.fa ${raw_reads_dr}/fq/${id}.fq.gz > mapping/${id}.longread.sam
fi

samtools sort -@ 39 mapping/${id}.longread.sam > mapping/${id}.longread.bam
bedtools coverage -sorted -a genome.overlappingwindows.bed -b mapping/${id}.longread.bam > mapping/genome.overlappingwindows.cov.tsv

# Gather metagenome pipeline results
echo "Gathering results in to results folder..."
# GC content
mv genome.GC.tsv results/
mv windows.GC.tsv results/
# rRNA predictions
mv genome.pro_blastn.rRNA.windows.bed results/
mv genome.euk_blastn.rRNA.windows.bed results/
# window blastn/x
mv mmseqs/${inDB}.${tag_pro_n}.bh results/
mv mmseqs/${inDB}.${tag_pro_x}.bh results/
mv mmseqs/${inDB}.${tag_euk_n}.bh results/
mv mmseqs/${inDB}.${tag_euk_x}.bh results/
mv mmseqs/${inDB}.${tag_euk_noAnts_n}.bh results/
mv mmseqs/${inDB}.${tag_human_n}.bh results/
mv mmseqs/*.m6 results/
mv genome.overlappingwindows.bed results/
mv genome.overlappingwindows.sorted.bed results/
mv mapping/genome.overlappingwindows.cov.tsv results/
# others
mv genome.overlappingwindows.tsv results
# clean up tmp files to save space
rm -rf tmp/*

# post processing results with a R script
cd results
Rscript /home/people/dinghe/github/GAGA-Metagenome-LGT/GAGA_metagenome_pipeline.R ${id}

###############################################################################################################
###  Lukas LGTfinder start
###############################################################################################################
# transform blast to bed
# m6 format: query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname
# -- blastn -- #
cat ${inDB}.${tag_pro_n}.m6|awk -F $'\t' 'BEGIN {OFS = FS} {if ($2<$3){print $1"\t"$2"\t"$3"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}else{print $1"\t"$3"\t"$2"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}}' > windows.fa.DB.pro.bed
cat ${inDB}.${tag_euk_n}.m6|awk -F $'\t' 'BEGIN {OFS = FS} {if ($2<$3){print $1"\t"$2"\t"$3"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}else{print $1"\t"$3"\t"$2"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}}' > windows.fa.DB.euk.bed

# Against how many species do we get hits?
cat windows.fa.DB.pro.bed|awk -F $'\t' 'BEGIN {OFS = FS} {if ($5>69) print $1,$2,$3,$4}'|perl -pe 's/\t.*\;/\t/g'|uniq |cut -f 1|uniq -c|perl -pe 's/ +([0-9]+) (.*)/$2\t$1/g'> hitcount.species.pro.tsv
cat windows.fa.DB.euk.bed|awk -F $'\t' 'BEGIN {OFS = FS} {if ($5>69) print $1,$2,$3,$4}'|perl -pe 's/\t.*\;/\t/g'|uniq |cut -f 1|uniq -c|perl -pe 's/ +([0-9]+) (.*)/$2\t$1/g'> hitcount.species.euk.tsv

#join files
LANG=en_EN join -i -1 1 -2 1 -a1 -a2 -o auto -e 0 -t $'\t'  <( LANG=en_EN sort -f -k1 hitcount.species.euk.tsv) <( LANG=en_EN sort -f -k1 hitcount.species.pro.tsv) > hitcount.species.tsv

# extract putative contaminants
cat hitcount.species.tsv|awk '{if ($2>0 && $2<=2 && $3>=100) print $0}' > putative.DB.contaminations.tsv

# -- blastx -- #
cat ${inDB}.${tag_pro_x}.m6|awk -F $'\t' 'BEGIN {OFS = FS} {if ($2<$3){print $1"\t"$2"\t"$3"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}else{print $1"\t"$3"\t"$2"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}}' > windows.fa.DB.prX.bed
cat ${inDB}.${tag_euk_x}.m6|awk -F $'\t' 'BEGIN {OFS = FS} {if ($2<$3){print $1"\t"$2"\t"$3"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}else{print $1"\t"$3"\t"$2"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}}' > windows.fa.DB.euX.bed

# -- Intersect euk and pro blast hits -- #
# get overlapping HSPs between pro (tmp1) and euk (tmp2).
# for each pro hit, all overlapping euk hits are returned (bedtools intersect -a windows.fa.DB.pro.bed -b windows.fa.DB.euk.bed -wao)
# sort overlap by euk bitscore (sort -t$'\t' -k 1,4 -k 10,10gr)
# keep only the best (-u) scoring euk HSP (as sorted before) (-k10,10gr) for each pro HSP ( sort -t$'\t' -u -k1,4)

# blastn
bedtools intersect -a windows.fa.DB.pro.bed -b windows.fa.DB.euk.bed -wao |sort -t$'\t' -k 1,4 -k 10,10gr | sort -t$'\t' -u -k1,4  > windows.fa.DB.proVSeuk.bed

# blastx
bedtools intersect -a windows.fa.DB.prX.bed -b windows.fa.DB.euX.bed -wao |sort -t$'\t' -k 1,4 -k 10,10gr | sort -t$'\t' -u -k1,4  > windows.fa.DB.proVSeuk.X.bed

# -- Filter 1 -- #
# Filter hsps where bacterial bitscore is higher than eukaryotic and where bacterial bitscore is > 69.
# bitscore 69 ~ evalue 1e-5
# bitscore 73 ~ evalue 1e-10

# blastn
cat windows.fa.DB.proVSeuk.bed |awk -F $'\t' '{if ($5>$10 && $5 > 69) print $0}'|sort -t$'\t' -k 1,3 -k 5,5gr | sort -t$'\t' -u -k1,3 > windows.fa.DB.LGTs.bed

# blastx
cat windows.fa.DB.proVSeuk.X.bed |awk -F $'\t' '{if ($5>$10 && $5 > 69) print $0}'|sort -t$'\t' -k 1,3 -k 5,5gr | sort -t$'\t' -u -k1,3 > windows.fa.DB.LGTs.X.bed

# -- Filter 2 -- #
# Filter LGT candidate windows
# keep only those where
# scaffold is not among the "contaminants"
# the bitscore difference is > 100
# the prokaryotic HSP is larger than 100 bp

## prepare grep-able list of contaminating scaffolds
### the output in the file "contaminants" needs to be all the contaminating scaffolds in the following patterns:
#### ^scaffold1017:
#### ^scaffold1161:
#### ^scaffold1288:
cat contaminantscaffolds.csv|sed 1d|cut -f1|perl -pe 's/^(.*)\n/^$1\:\n/g' > contaminants

#blastn
## strict filter: pro HSP length > 100 bp  & bitscore diff > 100
bedtools sort -i windows.fa.DB.LGTs.bed  |sort -k 1,4 -u|egrep -f contaminants - -v|awk -F $'\t' 'BEGIN {OFS = FS} {print $0"\t"$5-$10}'|sort -k 12,12gr|awk -F $'\t' 'BEGIN {OFS = FS} {if ($12 > 100 && sqrt(($2-$3)^2) > 100) print $0}' > windows.fa.DB.LGTs.filtered.bed
## loose filter: pro HSP length > 50 bp  & bitscore diff > 50
bedtools sort -i windows.fa.DB.LGTs.bed  |sort -k 1,4 -u|egrep -f contaminants - -v|awk -F $'\t' 'BEGIN {OFS = FS} {print $0"\t"$5-$10}'|sort -k 12,12gr|awk -F $'\t' 'BEGIN {OFS = FS} {if ($12 > 50 && sqrt(($2-$3)^2) > 50) print $0}' > windows.fa.DB.LGTs.filtered.loose.bed

#blastx
bedtools sort -i windows.fa.DB.LGTs.X.bed | sort -k 1,4 -u|egrep -f contaminants - -v | awk -F $'\t' 'BEGIN {OFS = FS} {print $0"\t"$5-$10}'|sort -k 12,12gr | awk -F $'\t' 'BEGIN {OFS = FS} {if ($12 > 100 && sqrt(($2-$3)^2) > 100) print $0}' > windows.fa.DB.LGTs.filtered.X.bed

# confirm contaminants
bedtools sort -i windows.fa.DB.LGTs.bed | sort -k 1,4 -u|egrep -f contaminants - | awk -F $'\t' 'BEGIN {OFS = FS} {print $0"\t"$5-$10}' | sort -k 12,12gr | awk -F $'\t' 'BEGIN {OFS = FS} {if ($12 < 0 ) print $0}' > windows.fa.DB.contaminants.filtered.bed


# -- merge HSPs to loci -- #
# split window coordinates in tabs (perl -pe 's/(^scaffold[0-9]+)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g')
# generate genomic position (awk  -F $'\t' 'BEGIN {OFS = FS} {print $1,$2+$4,$2+$5-1,$6,$7,$12,$14,$13}')
# merge loci <500 bp apart (bedtools merge -i - -c 4,5,6,7,8 -o first,collapse,collapse,collapse,collapse -d 500)
## report first pro hit (col 4)
## collapse bac bitscores (col 5)
## collapse euk bitscores (col 6)
## collapse pairwise differences in bitscores (col 7)
## collapse overlap (col 8)

# colnames of output
# scaffold  start end bestProHit;tstart;tend;evalue;bits;alnlen;pident;taxname; collapsed(proBit)  collapsed(eukBit)  collapsed(bitDiff) collapsed(overlapping)

# blastn
cat windows.fa.DB.LGTs.filtered.bed | perl -pe 's/(^.+?)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g'| awk  -F $'\t' 'BEGIN {OFS = FS} {print $1,$2+$4,$2+$5-1,$6,$7,$12,$14,$13}'| bedtools sort -i - | bedtools merge -i - -c 4,5,6,7,8 -o first,collapse,collapse,collapse,collapse -d 500 | bedtools sort -i - > LGTs.candidateloci.bed

cat windows.fa.DB.LGTs.filtered.loose.bed | perl -pe 's/(^.+?)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g'| awk  -F $'\t' 'BEGIN {OFS = FS} {print $1,$2+$4,$2+$5-1,$6,$7,$12,$14,$13}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6,7,8 -o first,collapse,collapse,collapse,collapse -d 100 | bedtools sort -i - > LGTs.candidateloci.loose.bed

# blastx
cat windows.fa.DB.LGTs.filtered.X.bed | perl -pe 's/(^.+?)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g' | awk  -F $'\t' 'BEGIN {OFS = FS} {print $1,$2+$4,$2+$5-1,$6,$7,$12,$14,$13}' | bedtools sort -i - |bedtools merge -i - -c 4,5,6,7,8 -o first,collapse,collapse,collapse,collapse -d 500 > LGTs.candidateloci.X.bed

# -- Analyses -- #
# extract regions to fasta files

bedtools getfasta -fi $base/genome.fa -bed LGTs.candidateloci.loose.bed -fo LGTs.candidateloci.loose.fa
bedtools getfasta -fi $base/genome.fa -bed LGTs.candidateloci.bed -fo LGTs.candidateloci.fa
bedtools getfasta -fi $base/genome.fa -bed LGTs.candidateloci.X.bed -fo LGTs.candidateloci.X.fa

# retrieve overlapping blastx hits
cat windows.fa.DB.prX.bed | perl -pe 's/(^.+?)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g' | awk -F $'\t' 'BEGIN {OFS = FS} {print $1,$2+$4,$2+$5-1,$6}' > loci.DB.prX.bed
cat loci.DB.prX.bed | sort -t$'\t' -u -k1,3 | bedtools sort -i -| bedtools merge -i - -c 4 -o collapse|bedtools intersect -a LGTs.candidateloci.bed -b - -wao|bedtools sort -i > LGTs.candidateloci.proteins.bed
cat loci.DB.prX.bed| sort -t$'\t' -u -k1,2 |bedtools sort -i -|bedtools merge -i - -c 4 -o collapse|bedtools intersect -a LGTs.candidateloci.loose.bed -b - -wao|bedtools sort -i > LGTs.candidateloci.loose.proteins.bed

## retrieve coverage in surrounding of loci
### plus/minus 20kb, i.e. 10 windows
cut -f 1,2 $base/genome.fa.fai > genome.file
bedtools slop -b 20000 -i LGTs.candidateloci.bed -g genome.file|bedtools intersect -b genome.overlappingwindows.cov.tsv -a - -wao |bedtools merge -i - -c 4,10,11,12,13,14 -o first,collapse,collapse,collapse,collapse,collapse|bedtools sort -i > LGTs.candidateloci.coverage.bed
bedtools slop -b 20000 -i LGTs.candidateloci.loose.bed -g genome.file|bedtools intersect -b genome.overlappingwindows.cov.tsv -a - -wao |bedtools merge -i - -c 4,10,11,12,13,14 -o first,collapse,collapse,collapse,collapse,collapse|bedtools sort -i > LGTs.candidateloci.loose.coverage.bed
bedtools slop -b 20000 -i LGTs.candidateloci.X.bed -g genome.file|bedtools intersect -b genome.overlappingwindows.cov.tsv -a - -wao |bedtools merge -i - -c 4,10,11,12,13,14 -o first,collapse,collapse,collapse,collapse,collapse|bedtools sort -i > LGTs.candidateloci.X.coverage.bed

# merge loci <5 kb apart
# blastn
bedtools merge -d 5000 -i LGTs.candidateloci.bed > LGTs.5kb.candidateregions.bed

# Retrieve subset from bam file for each candidate
## This could be improved to produce one bam file for each candidate region.
## Should be easy with parallel, but its too late to implement now.

bedtools intersect -abam $base/mapping/${id}.longread.bam -b LGTs.5kb.candidateregions.bed > LGTs.5kb.candidateregions.PacBio.bam
bedtools intersect -abam $base/mapping/${id}.longread.bam -b LGTs.candidateloci.loose.bed > LGTs.candidateloci.loose.PacBio.bam


## This should do the trick for all reads mapped to a contaminant scaffold
#cat contaminantscaffolds.csv|sed 1d|cut -f1|bedtools intersect -abam mapping/${id}.longread.bam -b - > contaminantScaffolds.PacBio.bam

#################################################################################################
# ## implement screen for low complexity
#################################################################################################

perl -I $toolsDir/SeqComplex/ $toolsDir/SeqComplex/profileComplexSeq.pl LGTs.candidateloci.loose.fa
perl -I $toolsDir/SeqComplex/ $toolsDir/SeqComplex/profileComplexSeq.pl LGTs.candidateloci.fa

###############################################################################################################
###  Lukas LGTfinder block ends
###############################################################################################################

# copy results to ERDA
tar -czvf ${inDB}.m6.tar.gz *.m6 --remove-files
lftp io.erda.dk -p 21 -e "mkdir -f /GAGA/Microbiome/Results/Latest/22012021/${id}; mirror -R $(pwd) /GAGA/Microbiome/Results/Latest/22012021/${id}; bye"

# Ending time/date
ENDTIME=$(date)
ENDTIME_INSEC=$(date +%s)
echo "==============================================="
echo "Pipeline started at $STARTTIME"
echo "Pipeline ended at $ENDTIME"
echo "Pipeline took $((ENDTIME_INSEC - STARTTIME_INSEC)) seconds to finish"
