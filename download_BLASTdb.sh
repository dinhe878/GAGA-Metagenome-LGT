#!/bin/bash

# load the compoterome modules

# go to the BLASTdb folder
cd /home/people/dinghe/ku_00039/people/dinghe/BLASTdb

# Generate a list of complete bacterial genome accession numbers from PATRIC database selection. This results in 1908 genomes and a total of 3471 accession numbers (chromosome & plasmid)
cat PATRIC_genome_list_21112020.txt | sed '1d' | awk -F '\t' '{print $0}' | cut -f 20 | sed 's/"//g' | sed 's/,/\n/g' > Prok_GenBankAcc.txt

# Split the list into 13 smaller list, each has 267 accession numbers
split -l 267 -d Prok_GenBankAcc.txt Prok_GenBankAcc.

# Generate qsub batch scripts for each list
for n in {00..12}
do

cat > retrieve_PATRIC_prok_geneomes.${n}.qsub <<EOF
#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00039 -A ku_00039
### Job name
#PBS -N PATRIC_prok_geneome.${n}
### Output files
#PBS -e PATRIC_prok_geneome.${n}.err
#PBS -o PATRIC_prok_geneome.${n}.log
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes/cores
#PBS -l nodes=1:ppn=1
### Memory
#PBS -l mem=10gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds>
#PBS -l walltime=1:12:00:00
### Set up the environmental variables
#PBS -V $NCBI_API_KEY=593466ec18bd4cdf7c4cfd2ea54ff3908209

# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

# Load all required modules for the job
module load tools perl edirect/7.50

### The real work begins here:
# Go to the working_dr
cd $HOME/ku_00039/people/dinghe/BLASTdb

# retrieve genomes
echo "dowdloading genomes..."
epost -db nuccore -input Prok_GenBankAcc.${n} -format acc | efetch -format fasta > PATRIC_1908_prok_geneomes.${n}.fa
echo "Genome dowdload completed"

# retrieve proteomes
echo "dowdloading proteomes..."
epost -db nuccore -input Prok_GenBankAcc.${n} -format acc | efetch -format fasta_cds_aa > PATRIC_1908_prok_geneomes_proteins.${n}.fa

echo "Proteome dowdload completed"

EOF

# Submit the jobs
echo "Submitting the job $PBS_JOBID for the list PATRIC_prok_geneome.${n}..."
qsub retrieve_PATRIC_prok_geneomes.${n}.qsub

done
