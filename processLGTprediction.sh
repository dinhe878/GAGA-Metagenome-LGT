#!/usr/bin/bash

####################################################################
# USAGE
####################################################################

usage()
{
    echo "usage: e.g. bash /Users/lukas/sciebo/Projects/LGT/GAGA-Metagenome-LGT/processLGTprediction.sh -s /Users/lukas/sciebo/Projects/LGT/GAGA-Metagenome-LGT/analyseLGTs.Rmd -b /Users/lukas/sciebo/Projects/LGT/results/ -i GAGA-0515 -d true"
    echo "-s: path to analyseLGTs.Rmd script."
    echo "-b: path to results folder."
    echo "-i: GAGA id to process."
    echo "-d: true=download from ftp; false=skip download."
}

####################################################################
# get options
####################################################################

while getopts s:b:i:d: flag
do
    case "${flag}" in
        s) export RMDpath=${OPTARG};;
        b) export base=${OPTARG};;
        i) export id=${OPTARG};;
        d) export download=${OPTARG};;
    esac
done


echo "script: $RMDpath";
echo "base: $base";
echo "id: ${id}";
echo "dowload: $download";

####################################################################
# Check if options are setup
####################################################################

if [[ -z "$RMDpath" || -z "$base" || -z "$id" || -z "$download"   ]]
then
   echo "Not all options set correctly.";
   usage;
   exit;
fi

####################################################################
#prepare directory
####################################################################

mkdir ${base}/${id}/
cd ${base}/${id}/

####################################################################
# Download data
####################################################################

{
if [ "$download" = true ]; then
# download data
lftp -p 21 io.erda.dk -e "\
set ftp:ssl-protect-data on;
set net:connection-limit 16;
lcd ${base}/${id}/;
cd ~/GAGA/Microbiome/Results/Latest/22012021/${id};
#basics
mirror --only-newer --file=genome.overlappingwindows.cov.tsv;
mirror --only-newer --file=genome.file;
mirror --only-newer --file=Taxa_screen.top5protaxa.pdf;
mirror --only-newer --file=/GAGA/Genome_assemblies/Final_PacBio_assemblies/${id}*.fasta;
# euk DB
mirror --only-newer --file=LGTs.candidateloci.loose.bed;
mirror --only-newer --file=LGTs.candidateloci.loose.coverage.bed;
mirror --only-newer --file=genome.overlappingwindows.cov.tsv;
mirror --only-newer --file=LGTs.candidateloci.loose.proteins.bed;
mirror --only-newer --file=LGTs.candidateloci.loose.fa;
mirror --only-newer --file=LGTs.candidateloci.loose.complex;
mirror --only-newer --file=LGTs.5kb.candidateregions.PacBio.bam;
mirror --only-newer --file=LGTs.candidateloci.loose.PacBio.bam;
mirror --only-newer --file=*.DB.euk_blastn.bh;
mirror --only-newer --file=*.DB.pro_blastn.bh;
# noAnts DB
mirror --only-newer --file=LGTs.nAo.candidateloci.loose.bed;
mirror --only-newer --file=LGTs.nAo.candidateloci.loose.coverage.bed;
mirror --only-newer --file=LGTs.nAo.candidateloci.loose.proteins.bed;
mirror --only-newer --file=LGTs.nAo.candidateloci.loose.fa;
mirror --only-newer --file=LGTs.nAo.candidateloci.loose.complex;
mirror --only-newer --file=LGTs.nAo.5kb.candidateregions.PacBio.bam;
mirror --only-newer --file=LGTs.nAo.candidateloci.loose.PacBio.bam;
mirror --only-newer --file=*.DB.euk_noAnts_blastn.bh;
bye;
"
fi
}

echo ""

############################################################
# check if all files are available
############################################################

{
if [ ! -f ${base}/${id}/LGTs.candidateloci.loose.bed ] || [ ! -f ${base}/${id}/LGTs.candidateloci.loose.coverage.bed ] || [ ! -f ${base}/${id}/genome.file ] ||  [ ! -f ${base}/${id}/${id}*.fasta ] ||  [ ! -f ${base}/${id}/LGTs.candidateloci.loose.PacBio.bam ]; then
    echo "Not all required files found!"
    echo "Files contained in ${base}/${id}/ are:"
    ls -lh ${base}/${id}/
    exit 0
fi
}

############################################################
# Run Rmarkdown script for this GAGA id
############################################################
# Rscript -e "rmarkdown::render('/Users/lukas/sciebo/Projects/LGT/GAGA-Metagenome-LGT/analyseLGTs.Rmd',output_file='/Users/lukas/sciebo/Projects/LGT/results/GAGA-0515.LGTfinder.html',params=list(id = 'GAGA-0515',dir='/Users/lukas/sciebo/Projects/LGT/results/',type='euk'))"
# Rscript -e "rmarkdown::render('/Users/lukas/sciebo/Projects/LGT/GAGA-Metagenome-LGT/analyseLGTs.Rmd',output_file='/Users/lukas/sciebo/Projects/LGT/results/GAGA-0515.LGTfinder.html',params=list(id = 'GAGA-0515',dir='/Users/lukas/sciebo/Projects/LGT/results/',type='noAnt'))"

Rscript -e "rmarkdown::render('${RMDpath}',output_file='${base}/${id}.LGTfinder.html',params=list(id = '${id}',dir='${base}',type='euk'))"
Rscript -e "rmarkdown::render('${RMDpath}',output_file='${base}/${id}.LGTfinder.html',params=list(id = '${id}',dir='${base}',type='noAnt'))"
