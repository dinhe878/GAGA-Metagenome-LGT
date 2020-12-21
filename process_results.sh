#!/bin/bash

#########################################################
# loading computerome modules
#########################################################

module load tools mmseqs2/release_12-113e3 emboss/6.6.0

#########################################################
# variables
#########################################################

results_dr=$HOME/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/GAGA-0200/results

#########################################################
# main processes
#########################################################

# processing prok results
# sort first by evalue, then by bitscore
cat $results_dr/prok/results_genome.m6 |sort -k1,1 -k7,7g -k8,8rg | perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > prok_genome.all.besthit.bls
cat $results_dr/prok/results_proteome.m6 |sort -k1,1 -k7,7g -k8,8rg | perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > prok_proteome.all.besthit.bls


# processing insect results
# sort first by evalue, then by bitscore
cat $results_dr/insect/resultsDB_insect.m6 |sort -k1,1 -k7,7g -k8,8rg | perl -pe 's/((.*?):(.*?)-(.*?)\t.*?)$/$2\t$3\t$4\t$1/g'|sort -k1,1 -k2,2g > insect.all.besthit.bls
