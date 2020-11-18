#!/bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -i GAGA-id -b parameterB -c parameterC"
   echo -e "\t-i GAGA-id"
   echo -e "\t-b Description of what is parameterB"
   echo -e "\t-c Description of what is parameterC"
   exit 1 # Exit script after printing help
}

while getopts "a:b:c:" opt
do
   case "$opt" in
      a ) GAGA-id="$OPTARG" ;;
      b ) parameterB="$OPTARG" ;;
      c ) parameterC="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -i "${GAGA-id}" ] || [ -z "$parameterB" ] || [ -z "$parameterC" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "${GAGA-id}"
echo "$parameterB"
echo "$parameterC"

#id=$1
#id=GAGA-0221
# set base directory for each genome to analyze
base=/home/people/dinghe/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${GAGA-id}}/
# get fasta.gz file for a given genome
# first batch of assemblies
genome=$(readlink -f /home/people/dinghe/ku_00039/people/dinghe/data/GAGA/Polished_assemblies/${GAGA-id}*_nextpolish.fasta.gz)
# second batch of assemblies (retrieved 12.06.2020)

# location of batch scripts
scripts=/home/people/dinghe/ku_00039/people/dinghe/scripts/batch/

mkdir $base
mkdir $base/genome
mkdir $base/genome/prok.bls
mkdir $base/genome/euk.bls
mkdir $base/genome/bac.rRNA
mkdir $base/genome/euk.rRNA
mkdir $base/err
mkdir $base/out
mkdir $base/results
mkdir $base/mmseqs
