1. Pull the recent git to CR2
```bash
cd /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/
git pull
```

2. Run the `LGTfinder.noAnts.sh` script on all genomes that have a `results` folder. This script compares blast results from the mmseq runs vs pro and the "noAnt" db to produce a set of LGT candidates. Candidates previously identified against the "euk" DB are ignored.
```bash
cd /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/

readlink -f /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/*/results| \
  perl -pe 's|.*\/(.*?)\/results|qsub -v \"id=$1"  /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTfinder.noAnts.sh|g' > /home/projects/ku_00039/people/luksch/GAGA/LGT/noAntRuns.qsub.sh

grep "GAGA-0275" /home/projects/ku_00039/people/luksch/GAGA/LGT/noAntRuns.qsub.sh > tmp.qsub
```

3. Check if the R scripts work on CR2.
Here, some R libraries need to be installed locally.
```bash
RMDpath=/home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/analyseLGTs.cr2.Rmd
outfolder=/home/projects/ku_00039/people/luksch/GAGA/LGT/
base=/home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/


cd $outfolder
module purge
module load tools gcc intel/perflibs pandoc/1.15 R/4.0.3
```

Install missing libraries.
```R
install.packages("openssl")
install.packages("randomcoloR")
BiocManager::install("rtracklayer")
BiocManager::install("bamsignals")
BiocManager::install("ggmsa")
```

4. Run the Rmarkdown script to filter and plot the LGT candidates.

```bash

readlink -f /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/*/results/LGTs.nAo.candidateloci.loose.fa| \
  perl -pe 's|.*\/(.*?)\/results/.*|qsub -v \"id=$1" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub|g' > /home/projects/ku_00039/people/luksch/GAGA/LGT/LGTplots.qsub.sh
bash /home/projects/ku_00039/people/luksch/GAGA/LGT/LGTplots.qsub.sh

#qsub -v "id=GAGA-0275" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
#qsub -v "id=GAGA-0001" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
#qsub -v "id=GAGA-0515" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=10:00:00" -v "id=GAGA-0002" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=10:00:00" -v "id=GAGA-0360" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=10:00:00" -v "id=GAGA-0396" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=5:00:00" -v "id=GAGA-0004" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=5:00:00" -v "id=GAGA-0014" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=5:00:00" -v "id=GAGA-0020" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=5:00:00" -v "id=GAGA-0024" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=5:00:00" -v "id=GAGA-0026" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=5:00:00" -v "id=GAGA-0028" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=5:00:00" -v "id=GAGA-0109" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=5:00:00" -v "id=GAGA-0114" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=5:00:00" -v "id=GAGA-0177" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=5:00:00" -v "id=GAGA-0187" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=5:00:00" -v "id=GAGA-0221" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=5:00:00" -v "id=GAGA-0346" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub
qsub -l "walltime=5:00:00" -v "id=GAGA-0454" /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTplots.qsub

```

The following PBS script is stored at `/home/projects/ku_00039/people/luksch/GAGA/LGT/LGTplots.qsub`
```bash
### Job name
#PBS -N LGTplots_${id}
### Output files
#PBS -e /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/run_log/LGTplots_${id}.err
#PBS -o /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/run_log/LGTplots_${id}.log
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes/cores
#PBS -l nodes=1:ppn=1:thinnode
### Minimum memory
#PBS -l mem=10gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds>
#PBS -l walltime=2:00:00

# qsub -v "id=GAGA-0275" /home/projects/ku_00039/people/luksch/GAGA/LGT/LGTplots.qsub

module load tools gcc intel/perflibs pandoc/1.15 R/4.0.3

#id=GAGA-0515
RMDpath=/home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/analyseLGTs.cr2.Rmd
outfolder=/home/projects/ku_00039/people/luksch/GAGA/LGT/
base=/home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/

Rscript -e "rmarkdown::render('${RMDpath}',output_file='${outfolder}/${id}.euk.LGTfinder.html',params=list(id = '${id}',dir='${base}',type='euk'))"
Rscript -e "rmarkdown::render('${RMDpath}',output_file='${outfolder}/${id}.noAnt.LGTfinder.html',params=list(id = '${id}',dir='${base}',type='noAnt'))"
```


Create tarball and download.
```bash

cd /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/



tar -cvzf /home/projects/ku_00039/people/luksch/GAGA/LGT/GAGA.LGT.tar.gz */results/LGTs* */results/lgt* */results/*.euk.lgt.* */results/*.noAnt.lgt.* */results/genome.file */results/Taxa* */results/contaminants*

scp luksch@ssh.computerome.dk:/home/projects/ku_00039/people/luksch/GAGA/LGT/GAGA-*.LGT.tar.gz /Users/lukas/sciebo/Projects/LGT/results/

```

Create tarball for each genome separately
```
qsub -v "id=GAGA-0001" /home/projects/ku_00039/people/luksch/GAGA/LGT/tarballing.qsub

readlink -f /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/*/results/*.noAnt.lgt.good.candidates.tsv| \
  perl -pe 's|.*\/(.*?)\/results/.*|qsub -v \"id=$1" /home/projects/ku_00039/people/luksch/GAGA/LGT/tarballing.qsub|g' > /home/projects/ku_00039/people/luksch/GAGA/LGT/tarballing.qsub.sh

bash /home/projects/ku_00039/people/luksch/GAGA/LGT/tarballing.qsub.sh
#scp luksch@ssh.computerome.dk:/home/projects/ku_00039/people/luksch/GAGA/LGT/GAGA-0515.LGT.tar.gz .

```

The following PBS script is stored at `/home/projects/ku_00039/people/luksch/GAGA/LGT/tarballing.qsub`
```bash
### Job name
#PBS -N tarball_${id}
### Output files
#PBS -e /home/projects/ku_00039/people/luksch/GAGA/LGT/log/tarball_${id}.err
#PBS -o /home/projects/ku_00039/people/luksch/GAGA/LGT/log/tarball_${id}.log
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes/cores
#PBS -l nodes=1:ppn=1:thinnode
### Minimum memory
#PBS -l mem=20gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds>
#PBS -l walltime=5:00:00

# qsub -v "id=GAGA-0275" /home/projects/ku_00039/people/luksch/GAGA/LGT/tarballing.qsub

tar -cvzf /home/projects/ku_00039/people/luksch/GAGA/LGT/${id}.LGT.tar.gz /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${id}/results/LGTs* /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${id}/results/lgt* /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${id}/results/*.euk.lgt.* /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${id}/results/*.noAnt.lgt.* /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${id}/results/genome.file /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${id}/results/Taxa* /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/${id}/results/contaminants*

```


Some genomes so far failed to run.
- GAGA-0002: R-script does not work properly for this genome.
- NCBI-0003,NCBI-0006,NCBI-0007: No euk LGT candidates detected
- NCBI-0010: no euk LGT candidates detected and job killed.
