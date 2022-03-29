# [Under construction!!]

# GAGA-Metagenome-LGT

This repository hosts a series of scripts to screen GAGA (Global Ant Genomics Alliance) assembled genomes for bacterial contigs/contaminants and LGT (Leteral Gene Transferred) gene candidates based on blast results (mmseqs blastn and blastx).

## Before running

The only tested working environment is on the HPC (High Performance Computing) clusters Computerome (computerome.dk) which is based on TORQUE queuing system (e.g. qsub) and Moab Workload Manager (e.g. msub). The dependencies (often requeired) will need to be loaded (e.g. module load; modules.sourceforge.net) within the HPC environment.

## The pipeline
In order to obtain the desired results, scripts need work as the following sequence.

1. runLGT.v4.sh
   * To run this script in the TORQUE queuing system with  necessary information:
   ```bash
   $ qsub -A <account_name> -l nodes=<node_count>:ppn=<core_count>:thinnode,mem=<mem_amount>,walltime=<time_requested> -N <Job Name> -e <your_run>.err -o <your_run>.log -v "id=<GAGA_id>,pacbio=<1 or ...>" runLGT.v4.sh
   ```
   * If want to receive email warning only when error occurs, add an option:
   ```bash
    -m n
   ```
2.

More are comming...
