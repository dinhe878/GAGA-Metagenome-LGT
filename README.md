# [Under construction!!]

# GAGA-Metagenome-LGT

This repository hosts a series of scripts to screen GAGA (Global Ant Genomics Alliance) assembled genomes for bacterial contigs/contaminants and LGT (lateral gene transfers) candidates based on blast results (mmseqs blastn and blastx).

## Before running

The only tested working environment is on the HPC (High Performance Computing) clusters Computerome (computerome.dk) which is based on TORQUE queuing system (e.g. qsub) and Moab Workload Manager (e.g. msub). The dependencies (often requeired) will need to be loaded (e.g. module load; modules.sourceforge.net) within the HPC environment.

## The pipeline

In order to obtain the desired results, scripts need work as the following sequence.

1. Run batch submission script [`runLGT.v4.sh`](qsub_batch/runLGT.v4.sh)
   * To run this script in the TORQUE queuing system with  necessary information:

   ```bash
   qsub -A <account_name> -l nodes=<node_count>:ppn=<core_count>:thinnode,mem=<mem_amount>,walltime=<time_requested> -N <Job Name> -e <your_run>.err -o <your_run>.log -v "id=<GAGA_id>,pacbio=<1 or ...>" runLGT.v4.sh
   ```

   * If want to receive email warning only when error occurs, add an option:

   ```bash
    -m n
   ```

   The script [`runLGT.v4.sh`](qsub_batch/runLGT.v4.sh) does the following for a given genome assembly:
   (A) blast sliding windows against different genome databases.  
   * generate 2.5 kb windows with 500 bp walking steps for all scaffolds
   * search each window against curated `mmseqs` genome databases, a *bacterial db*, an *insects db* with selected insect genome assemblies, an *no-ants insect db* with selected insect genome assemblies excluding ants, and a *human genome db*
   * keep the best hit for each window, sorted first by evalue (-k7,7g), then by bitscore (-k8,8gr)
   (B) Annotate rRNAs (with `barrnap`) in the assembly (to filter false-positive hits against the bacterial database)
   (C) Calculate gc-content and length (relevant only for windows at the end of a scaffold) for each window
   (D) Analyze sequencing raw read coverage (usually PacBio) for each window 
   
   </br>

2. The script invokes a separate R script ([`GAGA_metagenome_pipeline.R`](GAGA_metagenome_pipeline.R)) that processes the raw output (blast results, etc.) and identifies scaffolds to be flagged as bacterial or human contaminations.

   ```bash
   Rscript /home/people/dinghe/github/GAGA-Metagenome-LGT/GAGA_metagenome_pipeline.R ${id}
   ```

   Scaffolds are flagged according to the following criteria:

   ```
   # classify a scaffold as prokaryotic if
   #   1. 100 % of hitted windows are identified as prokaryotic OR
   #   2. more than 50% of windows are identified as prokaryotic AND bsratio > 200
   #   3. more than 50% of windows are identified as prokaryotic AND 0 < bsratio < 200 AND their
   #     3.1 bsratio.x > 50 AND
   #     3.2 GC content does not lie within 95% confidence interval of euk scaffolds GC content distribution AND
   #     3.3 coverage does not lie within 95% confidence interval of euk scaffolds coverage distribution
   # note that it is possible to have ratio 1 for only small amount of proWindow but the rest are uncertain (no hit)
   # for downstream genome annotation, recommending only the scaffolds tagged "pro" should be excluded
   # for downstram microbiome analyses, recommending to include all "pro" and "uncertain" scaffolds
   # Checking for human contamination has added: kingdom will be assigned human if more than half of
   # the windows hitted better against human DB
   ```

3. To identify regions in the genomes that have the signature of being horizontally transferred from a bacterium, the script [`runLGT.v4.sh`](qsub_batch/runLGT.v4.sh) does the following:
   ```
   # -- Intersect euk and pro blast hits -- #
   # get overlapping HSPs between pro (windows.fa.DB.pro.bed) and euk (windows.fa.DB.euk.bed).
   # for each pro hit, all overlapping euk hits are returned (bedtools intersect -a windows.fa.DB.pro.bed -b windows.fa.DB.euk.bed -wao)
   # sort overlap by euk bitscore (sort -t$'\t' -k 1,4 -k 10,10gr)
   # keep only the best (-u) scoring euk HSP (as sorted before) (-k10,10gr) for each pro HSP ( sort -t$'\t' -u -k1,4)
   ```

More are comming...
