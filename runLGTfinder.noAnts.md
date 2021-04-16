```bash
cd /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/
git pull

cd /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/

readlink -f /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/*/results| \
  perl -pe 's|.*\/(.*?)\/results|qsub -v \"id=$1"  /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTfinder.noAnts.sh|g' > /home/projects/ku_00039/people/luksch/GAGA/LGT/noAntRuns.qsub.sh

grep "GAGA-0275" /home/projects/ku_00039/people/luksch/GAGA/LGT/noAntRuns.qsub.sh > tmp.qsub

RMDpath=/home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/analyseLGTs.cr2.Rmd
outfolder=/home/projects/ku_00039/people/luksch/GAGA/LGT/
base=/home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/


cd $outfolder
module purge
module load tools gcc intel/perflibs pandoc/1.15 R/4.0.3
```


```R
install.packages("openssl")
install.packages("randomcoloR")
BiocManager::install("rtracklayer")
BiocManager::install("bamsignals")
BiocManager::install("ggmsa")
```

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

# qsub -v "id=GAGA-0515" /home/projects/ku_00039/people/luksch/GAGA/LGT/LGTplots.qsub

module load tools gcc intel/perflibs pandoc/1.15 R/4.0.3

#id=GAGA-0515
RMDpath=/home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/analyseLGTs.cr2.Rmd
outfolder=/home/projects/ku_00039/people/luksch/GAGA/LGT/
base=/home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/

Rscript -e "rmarkdown::render('${RMDpath}',output_file='${outfolder}/${id}.euk.LGTfinder.html',params=list(id = '${id}',dir='${base}',type='euk'))"
Rscript -e "rmarkdown::render('${RMDpath}',output_file='${outfolder}/${id}.noAnt.LGTfinder.html',params=list(id = '${id}',dir='${base}',type='noAnt'))"
```



```bash

cd ${base}
tar -cvzf $outfolder/GAGA.LGT.tar.gz */results/LGTs* */results/lgt* */results/${id}.euk.lgt.* */results/${id}.noAnt.lgt.* */results/genome.file

scp luksch@ssh.computerome.dk:/home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/*/results/*LGT.tar.gz /Users/lukas/sciebo/Projects/LGT/results/

```
