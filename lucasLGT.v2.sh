#!/bin/bash

#########################################################
# loading computerome modules
#########################################################

module load tools samtools/1.10 bedtools/2.28.0 pigz/2.3.4 seqkit/0.13.2

#########################################################
# commandline input options
#########################################################

helpFunction()
{
   echo ""
   echo "Usage: $0 -i GAGA-id"
   echo -e "\t-i GAGA-id"
   exit 1 # Exit script after printing help
}

while getopts "i:" opt
do
   case "$opt" in
      i ) id="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$id" ]
then
   echo "Parameters are missing";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "Now processing $id..."

#########################################################
# setup variables and folder structure
#########################################################

# set base directory for each genome to analyze
base=/home/people/dinghe/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/$id/

# get fasta.gz file for a given genome
# Polished assemblies
genome=$(readlink -f /home/people/dinghe/ku_00039/people/dinghe/data/GAGA/Polished_assemblies/$id*_nextpolish.fasta.gz)

# location of batch scripts
scripts=/home/people/dinghe/ku_00039/people/dinghe/scripts/batch/

mkdir $base
cd $base
mkdir $base/genome
mkdir $base/genome/prok.bls
mkdir $base/genome/euk.bls
mkdir $base/genome/bac.rRNA
mkdir $base/genome/euk.rRNA
mkdir $base/err
mkdir $base/out
mkdir $base/results
mkdir $base/mmseqs

#########################################################
# prepare input files
#########################################################

# prepare genome fasta
zcat $genome|cut -f 1 -d "|"  > $base/genome.fa
samtools faidx genome.fa

# generate 2.5 kb windows with 500 bp walking steps (overlap 500 bp with the previous and the subsequent window)
seqkit sliding -s 2000 -W 2500 -g genome.fa > windows.fa
