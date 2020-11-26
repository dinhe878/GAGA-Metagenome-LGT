#!/bin/bash

# go to the BLASTdb folder
cd /home/people/dinghe/ku_00039/people/dinghe/BLASTdb/EDirect

# make directories
mkdir Prok_GenBankAcc
acc_dir=./Prok_GenBankAcc
mkdir Prok_Genome_fasta
genome_dir=./Prok_Genome_fasta
mkdir Prok_Proteome_fasta
proteome_dir=./Prok_Proteome_fasta

# Generate a list of complete bacterial genome accession numbers from PATRIC database selection. This results in 1908 genomes and a total of 3471 accession numbers (chromosome & plasmid)
cat PATRIC_genome_list_21112020.txt | sed '1d' | awk -F '\t' '{print $0}' | cut -f 20 | sed 's/"//g' | sed 's/,/\n/g' > Prok_GenBankAcc.txt

# Split the list into smaller lists of 10 accession numbers due to api_key access limit
split -l 10 -d Prok_GenBankAcc.txt $acc_dir/Prok_GenBankAcc.

# Generate qsub a batch script
cat > retrieve_PATRIC_prok_geneomes.qsub <<EOF
#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00039 -A ku_00039
### Job name
#PBS -N PATRIC_prok_geneome
### Output files
#PBS -e PATRIC_prok_geneome.err
#PBS -o PATRIC_prok_geneome.log
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes/cores
#PBS -l nodes=1:ppn=1
### Memory
#PBS -l mem=10gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds>
#PBS -l walltime=10:00:00
### Set up the environmental variables
NCBI_API_KEY=593466ec18bd4cdf7c4cfd2ea54ff3908209

# Go to the directory from where the job was submitted (initial directory is \$HOME)
echo "Working directory is \$PBS_O_WORKDIR"
cd \$PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS=\`wc -l < \$PBS_NODEFILE\`
echo This job has allocated \$NPROCS nodes

# Load all required modules for the job
module load tools perl edirect/7.50

### The real work begins here:
# Go to the working_dr
cd \$HOME/ku_00039/people/dinghe/BLASTdb/EDirect
acc_dir=./Prok_GenBankAcc

# retrieve genomes
for f in \$acc_dir
do
  echo "dowdloading genomes..."
  epost -db nuccore -input \$f -format acc -api_key \$NCBI_API_KEY | efetch -format fasta > \$f.fa
  echo \$f

  # try to avoid spamming the NCBI server
  sleep 2

  # retrieve proteomes
  echo "dowdloading proteomes..."
  epost -db nuccore -input \$f -format acc -api_key \$NCBI_API_KEY | efetch -format fasta_cds_aa > \$f.proteins.fa
  echo \$f
done

EOF

# Submit the jobs
qsub retrieve_PATRIC_prok_geneomes.qsub
echo "Batch job submitted"
