#!/usr/bin/bash

## usage:
## bash /Users/lukas/sciebo/Projects/LGT/GAGA-Metagenome-LGT/runLGTfinder.sh /Users/lukas/sciebo/Projects/LGT/GAGA-Metagenome-LGT/analyseLGTs.Rmd /Users/lukas/sciebo/Projects/LGT/results/ GAGA-0024
#Define base dir and target Rmarkdown script
#RMDpath=/Users/lukas/sciebo/Projects/LGT/GAGA-Metagenome-LGT/analyseLGTs.Rmd
#base=/Users/lukas/sciebo/Projects/LGT/results/

# Define GAGA id to analyze
#export id=GAGA-0003

if [ $# -ne 3 ]
  then
    echo "
    Not all arguments supplied.

    usage:
    e.g. bash /Users/lukas/sciebo/Projects/LGT/GAGA-Metagenome-LGT/processLGTprediction.sh /Users/lukas/sciebo/Projects/LGT/GAGA-Metagenome-LGT/analyseLGTs.Rmd /Users/lukas/sciebo/Projects/LGT/results/ GAGA-0024
    "
    exit 0
fi


# to run from command line
RMDpath=$1
base=$2
id=$3

#prepare directory
mkdir ${base}/${id}/
cd ${base}/${id}/
{
if [ ! -f ${base}/${id}/LGTs.candidateloci.loose.bed ]; then
# download data
lftp -p 21 io.erda.dk -e "\
set ftp:ssl-protect-data on;
set net:connection-limit 16;
lcd ${base}/${id}/
mget GAGA/Microbiome/Results/Latest/22012021/${id}/LGTs.candidateloci.*;
mget GAGA/Microbiome/Results/Latest/22012021/${id}/LGTs.5kb.*;
mget GAGA/Microbiome/Results/Latest/22012021/${id}/genome.overlappingwindows.cov.tsv;
mget GAGA/Microbiome/Results/Latest/22012021/${id}/*.pdf;
mget GAGA/Microbiome/Results/Latest/22012021/${id}/*.tsv;
mget GAGA/Genome_assemblies/Final_PacBio_assemblies/${id}*.fasta;
bye;
"
fi
}

echo ""
############################################################
# check if all files are available
############################################################

{
if [ ! -f ${base}/${id}/LGTs.candidateloci.loose.bed ]; then
    echo "File LGTs.candidateloci.loose.bed not found!"
    echo "Files contained in ${base}/${id}/ are:"
    ls -lh ${base}/${id}/
    exit 0
fi
}

{
if [ ! -f ${base}/${id}/LGTs.candidateloci.loose.coverage.bed ]; then
    echo "File LGTs.candidateloci.loose.coverage.bed not found!"
    echo "Files contained in ${base}/${id}/ are:"
    ls -lh ${base}/${id}/
    exit 0
fi
}

{
if [ ! -f ${base}/${id}/genome.overlappingwindows.cov.tsv ]; then
    echo "File genome.overlappingwindows.cov.tsv not found!"
    echo "Files contained in ${base}/${id}/ are:"
    ls -lh ${base}/${id}/
    exit 0
fi
}

{
if [ ! -f ${base}/${id}/LGTs.candidateloci.loose.proteins.bed ]; then
    echo "File LGTs.candidateloci.loose.proteins.bed not found!"
    echo "Files contained in ${base}/${id}/ are:"
    ls -lh ${base}/${id}/
    exit 0
fi
}

{
if [ ! -f ${base}/${id}/LGTs.candidateloci.loose.fa ]; then
    echo "File LGTs.candidateloci.loose.fa not found!"
    echo "Files contained in ${base}/${id}/ are:"
    ls -lh ${base}/${id}/
    exit 0
fi
}

{
if [ ! -f ${base}/${id}/LGTs.candidateloci.loose.complex ]; then
    echo "File LGTs.candidateloci.loose.complex not found!"
    echo "Files contained in ${base}/${id}/ are:"
    ls -lh ${base}/${id}/
    exit 0
fi
}

{
if [ ! -f ${base}/${id}/LGTs.5kb.candidateregions.PacBio.bam ]; then
    echo "File LGTs.5kb.candidateregions.PacBio.bam not found!"
    echo "Files contained in ${base}/${id}/ are:"
    ls -lh ${base}/${id}/
    exit 0
fi
}

############################################################
# Run Rmarkdown script for this GAGA id
############################################################

Rscript -e "rmarkdown::render('${RMDpath}',output_file='${base}/${id}.LGTfinder.html',params=list(id = '${id}',dir='${base}'))"
