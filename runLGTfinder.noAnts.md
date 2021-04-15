```bash
cd /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/
git pull

cd /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/

readlink -f /home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/*/results| \
  perl -pe 's|.*\/(.*?)\/results|qsub -v \"id=$1"  /home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/LGTfinder.noAnts.sh|g' > ~/GAGA/LGT/noAntRuns.qsub.sh

RMDpath=/home/projects/ku_00039/people/luksch/software/GAGA-Metagenome-LGT/analyseLGTs.cr2.Rmd
outfolder=/home/projects/ku_00039/people/luksch/GAGA/LGT/
base=/home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/
id=GAGA-0001

cd $outfolder
module purge
module load tools gcc intel/perflibs
module load pandoc/1.15 R/4.0.3
```


```R
install.packages("openssl")
install.packages("randomcoloR")
BiocManager::install("rtracklayer")
BiocManager::install("bamsignals")
BiocManager::install("ggmsa")
```

```bash
Rscript -e "rmarkdown::render('${RMDpath}',output_file='${outfolder}/${id}.LGTfinder.html',params=list(id = '${id}',dir='${base}',type='euk'))"
Rscript -e "rmarkdown::render('${RMDpath}',output_file='${outfolder}/${id}.LGTfinder.html',params=list(id = '${id}',dir='${base}',type='noAnt'))"
```
