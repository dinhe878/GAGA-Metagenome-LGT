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

module load tools perl samtools/1.10 bedtools/2.28.0 pigz/2.3.4 mmseqs2/release_12-113e3 barrnap/0.7 emboss/6.6.0 minimap2/2.17r941 gcc intel/perflibs R/3.6.1 lftp/4.9.2

#########################################################
# setup variables and folder structure                  #
#########################################################

# Starting time/date
STARTTIME=$(date)
STARTTIME_INSEC=$(date +%s)

# GAGA-ID is passed through commandline option, e.g. -v "id=GAGA-0024"

# set base directory for each genome to analyze
base=/home/people/dinghe/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${id}/

# set variables pointing to the final GAGA genome assembly
assembly_dr=/home/people/dinghe/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Final_PacBio_assemblies/
genome_file_name=$(ls -l $assembly_dr | awk -v pat="${id}" '$0~pat' | awk '{split($0,a," "); print a[9]}')
genome=${assembly_dr}${genome_file_name}

# GAGA genome pacbio raw reads folder
raw_reads_dr=/home/people/dinghe/ku_00039/people/dinghe/data/GAGA/Raw_genome_reads

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
cat genome.overlappingwindows.bed|perl -pe 's/(.*?)\t(.*?)\t(.*?)/$1\:$2\-$3/g'  > genome.overlappingwindows.tsv
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
genome_db_euk="$targetBlastnDB/mmseq.genome.71_clean_insectDB"
proteome_db_euk="$targetBlastxDB/Insecta"

# run mmseqs search twice, first on prokaryotic (i.e. tag="pro") and then on eukaryotic (tag="euk")
echo "Starting mmseqs blasting ${inDB} against ${genome_db_pro}..."
mmseqs search ${inDB} ${genome_db_pro} mmseqs/${inDB}.${tag_pro_n}.resDB tmp --start-sens 1 --sens-steps 2 -s ${sensitivity} --search-type 3 1> mmseqs/${inDB}.${tag_pro_n}.search.out 2> mmseqs/${inDB}.${tag_pro_n}.search.err
echo "Starting mmseqs blasting ${inDB} against ${proteome_db_pro}..."
mmseqs search ${inDB} ${proteome_db_pro} mmseqs/${inDB}.${tag_pro_x}.resDB tmp --start-sens 1 --sens-steps 2 -s ${sensitivity} --search-type 3 1> mmseqs/${inDB}.${tag_pro_x}.search.out 2> mmseqs/${inDB}.${tag_pro_x}.search.err
echo "Starting mmseqs blasting ${inDB} against ${genome_db_euk}..."
mmseqs search ${inDB} ${genome_db_euk} mmseqs/${inDB}.${tag_euk_n}.resDB tmp --start-sens 1 --sens-steps 2 -s ${sensitivity} --search-type 3 1> mmseqs/${inDB}.${tag_euk_n}.search.out 2> mmseqs/${inDB}.${tag_euk_n}.search.err
echo "Starting mmseqs blasting ${inDB} against ${proteome_db_euk}..."
mmseqs search ${inDB} ${proteome_db_euk} mmseqs/${inDB}.${tag_euk_x}.resDB tmp --start-sens 1 --sens-steps 2 -s ${sensitivity} --search-type 3 1> mmseqs/${inDB}.${tag_euk_x}.search.out 2> mmseqs/${inDB}.${tag_euk_x}.search.err

# run mmseqs convertalis (to generate blast m6-like output format)
echo "Converting blast results..."
mmseqs convertalis ${inDB} ${genome_db_pro} mmseqs/${inDB}.${tag_pro_n}.resDB mmseqs/${inDB}.${tag_pro_n}.m6 --format-output query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname  1> mmseqs/${inDB}.${tag_pro_n}.convert.out 2> mmseqs/${inDB}.${tag_pro_n}.convert.err
mmseqs convertalis ${inDB} ${proteome_db_pro} mmseqs/${inDB}.${tag_pro_x}.resDB mmseqs/${inDB}.${tag_pro_x}.m6 --format-output query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname  1> mmseqs/${inDB}.${tag_pro_x}.convert.out 2> mmseqs/${inDB}.${tag_pro_x}.convert.err

mmseqs convertalis ${inDB} ${genome_db_euk} mmseqs/${inDB}.${tag_euk_n}.resDB mmseqs/${inDB}.${tag_euk_n}.m6 --format-output query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname  1> mmseqs/${inDB}.${tag_euk_n}.convert.out 2> mmseqs/${inDB}.${tag_euk_n}.convert.err
mmseqs convertalis ${inDB} ${proteome_db_euk} mmseqs/${inDB}.${tag_euk_x}.resDB mmseqs/${inDB}.${tag_euk_x}.m6 --format-output query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname  1> mmseqs/${inDB}.${tag_euk_x}.convert.out 2> mmseqs/${inDB}.${tag_euk_x}.convert.err

# sort first by evalue (-k7,7g), then by bitscore (-k8,8gr)
# keep best hit only
cat mmseqs/${inDB}.${tag_pro_n}.m6 |sort -k1,1 -k7,7g -k8,8gr | sort -u -k1,1 --merge |perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > mmseqs/${inDB}.${tag_pro_n}.bh
cat mmseqs/${inDB}.${tag_pro_x}.m6 |sort -k1,1 -k7,7g -k8,8gr | sort -u -k1,1 --merge |perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > mmseqs/${inDB}.${tag_pro_x}.bh
cat mmseqs/${inDB}.${tag_euk_n}.m6 |sort -k1,1 -k7,7g -k8,8gr | sort -u -k1,1 --merge |perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > mmseqs/${inDB}.${tag_euk_n}.bh
cat mmseqs/${inDB}.${tag_euk_x}.m6 |sort -k1,1 -k7,7g -k8,8gr | sort -u -k1,1 --merge |perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > mmseqs/${inDB}.${tag_euk_x}.bh

# retrieve a list of all the windows that hit against the prokDB
cat mmseqs/${inDB}.${tag_pro_n}.bh |cut -f 4 > mmseqs/${inDB}.${tag_pro_n}.bh.lst
cat mmseqs/${inDB}.${tag_pro_x}.bh |cut -f 4 > mmseqs/${inDB}.${tag_pro_x}.bh.lst
cat mmseqs/${inDB}.${tag_euk_n}.bh |cut -f 4 > mmseqs/${inDB}.${tag_euk_n}.bh.lst
cat mmseqs/${inDB}.${tag_euk_x}.bh |cut -f 4 > mmseqs/${inDB}.${tag_euk_x}.bh.lst

# rRNA prediction
echo "Starting rRNA prediction..."
barrnap --threads 40 -k bac genome.fa > genome.${tag_pro_n}.rRNA.gff3
barrnap --threads 40 -k euk genome.fa > genome.${tag_euk_n}.rRNA.gff3
## retrieve windows that overlap an rRNA
bedtools intersect -a genome.overlappingwindows.bed -b genome.${tag_pro_n}.rRNA.gff3 -wa -wb > genome.${tag_pro_n}.rRNA.windows.bed
bedtools intersect -a genome.overlappingwindows.bed -b genome.${tag_euk_n}.rRNA.gff3 -wa -wb > genome.${tag_euk_n}.rRNA.windows.bed

# Calculate GC content and length for each scaffold
echo "Calculating GC content of scaffolds..."
infoseq  -nocolumn -delimiter "\t" -auto -only -name -length -pgc genome.fa > genome.GC.tsv

# calculate gc and length by window
cat windows.fa|perl -pe 's/(>.*?)\:(.*)/$1\@$2/g' > tmp.fa
infoseq  -nocolumn -delimiter "\t" -auto -only -name -length -pgc tmp.fa |perl -pe 's/^(.*?)@(.*)/$1:$2/g' > ./windows.GC.tsv
rm tmp.fa

# Analyze long-read coverage
echo "Gathering sequencing coverage information..."
bamToFastq -i ${raw_reads_dr}/bam/${id}.bam -fq ${raw_reads_dr}/fq/${id}.fq.gz
minimap2 -t 40 -ax map-pb genome.fa ${raw_reads_dr}/fq/${id}.fq.gz > mapping/${id}.longread.sam
samtools view -S -b mapping/${id}.longread.sam | samtools sort > mapping/${id}.longread.bam
bedtools coverage -a genome.overlappingwindows.bed -b mapping/${id}.longread.bam > mapping/genome.overlappingwindows.cov.tsv
rm ${raw_reads_dr}/fq/${id}.fq.gz

# Gather results
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
mv mmseqs/*.m6 results/
mv genome.overlappingwindows.bed results/
mv mapping/genome.overlappingwindows.cov.tsv results/
# others
mv genome.overlappingwindows.tsv results
# clean up tmp files to save space
rm -rf tmp/*

# post processing results with a R script
cd results
Rscript /home/people/dinghe/github/GAGA-Metagenome-LGT/GAGA_metagenome_pipeline.R ${id}

# copy results to ERDA
tar -czvf ${inDB}.m6.tar.gz *.m6
lftp io.erda.dk -p 21 -e "mkdir -f /GAGA/Microbiome/Results/Latest/22012021/${id}; mirror -R $(pwd) /GAGA/Microbiome/Results/Latest/22012021/${id}; bye"

# Ending time/date
ENDTIME=$(date)
ENDTIME_INSEC=$(date +%s)
echo "==============================================="
echo "Pipeline started at $STARTTIME"
echo "Pipeline ended at $ENDTIME"
echo "Pipeline took $((ENDTIME_INSEC - STARTTIME_INSEC)) seconds to finish"
