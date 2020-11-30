#!/bin/bash

# go to the BLASTdb folder
cd $HOME/ku_00039/people/dinghe/BLASTdb/EDirect

# make directories
mkdir Prok_GenBankAcc
acc_dir=./Prok_GenBankAcc
mkdir Prok_Genome_fasta
prok_genome_dir=./Prok_Genome_fasta
mkdir Prok_Proteome_fasta
prok_proteome_dir=./Prok_Proteome_fasta
mkdir Insect_Genome_fasta
insect_genome_dir=./Insect_Genome_fasta
mkdir Insect_Proteome_fasta
insect_proteome_dir=./Prok_Proteome_fasta


# Generate a list of complete bacterial genome accession numbers from PATRIC database selection. This results in 1908 genomes and a total of 3471 accession numbers (chromosome & plasmid)
cat PATRIC_genome_list_21112020.txt | sed '1d' | awk -F '\t' '{print $0}' | cut -f 20 | sed 's/"//g' | sed 's/,/\n/g' > Prok_GenBankAcc.txt
# Split the list into smaller lists of 10 accession numbers due to api_key access limit
split -l 10 -d Prok_GenBankAcc.txt $acc_dir/Prok_GenBankAcc.

# Generate a list of complete insect genome/chromosome accession numbers from NCBI database. This results in 43 genomes.
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt
cat eukaryotes.txt | awk -F '\t' '{if ($17=="Chromosome" || $17=="Complete Genome") print $0}'| awk -F ' ' '!seen[$1]++'|awk -F '\t' '{if ($6=="Insects") print $0}' | cut -f 9 | sed 's/"//g' | sed 's/,/\n/g' > Insect_GenBankAcc.txt


# Generate qsub a batch script
cat > retrieve_geneomes.qsub <<EOF
#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00039 -A ku_00039
### Job name
#PBS -N PATRIC_prok_geneome
### Output files
#PBS -e Download_geneome.err
#PBS -o Download_geneome.log
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
# Go to the prok_working_dr
cd \$HOME/ku_00039/people/dinghe/BLASTdb/EDirect
acc_dir=\$HOME/ku_00039/people/dinghe/BLASTdb/EDirect/Prok_GenBankAcc
prok_genome_dir=\$HOME/ku_00039/people/dinghe/BLASTdb/EDirect/Prok_Genome_fasta
prok_proteome_dir=\$HOME/ku_00039/people/dinghe/BLASTdb/EDirect/Prok_Proteome_fasta
insect_genome_dir=\$HOME/ku_00039/people/dinghe/BLASTdb/EDirect/Insect_Genome_fasta
insect_proteome_dir=\$HOME/ku_00039/people/dinghe/BLASTdb/EDirect/Insect_Proteome_fasta

# retrieve prok genomes/proteomes
while:
epost -db nuccore -input Insect_GenBankAcc.txt -format acc -api_key \$NCBI_API_KEY | efetch -format fasta > Insect_34_genome.fa || { break }
sleep 2
True
done

while:
epost -db nuccore -input Insect_GenBankAcc.txt -format acc -api_key \$NCBI_API_KEY | efetch -format fasta_cds_aa > Insect_34_proteins.fa || { break }
sleep 2
True
done

# retrieve prok genomes/proteomes
for f in \$acc_dir/*
do
  fName=\`basename \$f\`
  echo \$fName

  echo "dowdloading prok genomes..."
  while:
  epost -db nuccore -input \$acc_dir/\$fName -format acc -api_key \$NCBI_API_KEY | efetch -format fasta > \$fName.fa || { break }
  True
  done
  mv \$fName.fa \$prok_genome_dir

  # try to avoid spamming the NCBI server
  sleep 2

  echo "dowdloading prok proteomes..."
  while:
  epost -db nuccore -input \$acc_dir/\$fName -format acc -api_key \$NCBI_API_KEY | efetch -format fasta_cds_aa > \$fName.proteins.fa || { break }
  True
  done
  mv \$fName.proteins.fa \$prok_proteome_dir

done

EOF

# Submit the jobs
#qsub retrieve_geneomes.qsub
#echo "Batch job submitted"
