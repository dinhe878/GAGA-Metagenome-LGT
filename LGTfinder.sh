#cd /global/scratch/schradel/BLAST/LGTv3/


#################################################################################################
# Define query
#################################################################################################
#query=GAGA-0001
query=GAGA-0024
query=GAGA-0221
id=GAGA-0024
#################################################################################################

#################################################################################################
# switch to results folder
#################################################################################################
cd /global/scratch/schradel/BLAST/LGTv3/$query/
base=/home/people/dinghe/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA/$id/
cd $base
#################################################################################################

#################################################################################################
# transform blast to bed
#################################################################################################
### m6 format: query,qstart,qend,target,tstart,tend,evalue,bits,alnlen,pident,taxlineage,taxid,taxname

## blastn
cat query.$query.DB.pro_blastn.m6|awk -F $'\t' 'BEGIN {OFS = FS} {if ($2<$3){print $1"\t"$2"\t"$3"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}else{print $1"\t"$3"\t"$2"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}}' > windows.fa.DB.pro.bed
cat query.$query.DB.euk_blastn.m6|awk -F $'\t' 'BEGIN {OFS = FS} {if ($2<$3){print $1"\t"$2"\t"$3"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}else{print $1"\t"$3"\t"$2"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}}' > windows.fa.DB.euk.bed

# Against how many species do we get hits?
cat windows.fa.DB.pro.bed|awk -F $'\t' 'BEGIN {OFS = FS} {if ($5>69) print $1,$2,$3,$4}'|perl -pe 's/\t.*\;/\t/g'|uniq |cut -f 1|uniq -c|perl -pe 's/ +([0-9]+) (.*)/$2\t$1/g'> hitcount.species.pro.tsv
cat windows.fa.DB.euk.bed|awk -F $'\t' 'BEGIN {OFS = FS} {if ($5>69) print $1,$2,$3,$4}'|perl -pe 's/\t.*\;/\t/g'|uniq |cut -f 1|uniq -c|perl -pe 's/ +([0-9]+) (.*)/$2\t$1/g'> hitcount.species.euk.tsv

#join files
LANG=en_EN join -i -1 1 -2 1 -a1 -a2 -o auto -e 0 -t $'\t'  <( LANG=en_EN sort -f -k1 hitcount.species.euk.tsv) <( LANG=en_EN sort -f -k1 hitcount.species.pro.tsv) > hitcount.species.tsv

# extract putative contaminants
cat hitcount.species.tsv|awk '{if ($2>0 && $2<=2 && $3>=100) print $0}' > putative.DB.contaminations.tsv


#cat windows.fa.DB.euk.bed|awk '{if ($5>69) print $1}'|uniq -c|perl -pe 's/ +([0-9]+) (.*)/$2\t$1/g'> hitcount.euk.tsv
#cat windows.fa.DB.pro.bed|awk '{if ($5>69) print $1}'|uniq -c|perl -pe 's/ +([0-9]+) (.*)/$2\t$1/g'> hitcount.pro.tsv



grep "scaffold29:1777500-1779500" windows.fa.DB.euk.bed

bedtools merge -i windows.fa.DB.pro.bed -c 5 -o stdev

# quantify drop in evalue
 cat windows.fa.DB.pro.bed|perl -pe 's/(^.+?)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g'|bedtools sort -i -|bedtools merge -i - -c 5 -o sstdev|head

cat hitcount.euk.tsv|sort > tmp1
cat hitcount.pro.tsv|sort > tmp2

# eliminate things that only hit against one to three references
# eliminate things that only have ONE very good hit
# quantify the drop in evalue
## blastx
cat query.$query.DB.pro_blastx.m6|awk -F $'\t' 'BEGIN {OFS = FS} {if ($2<$3){print $1"\t"$2"\t"$3"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}else{print $1"\t"$3"\t"$2"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}}' > windows.fa.DB.prX.bed
cat query.$query.DB.euk_blastx.m6|awk -F $'\t' 'BEGIN {OFS = FS} {if ($2<$3){print $1"\t"$2"\t"$3"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}else{print $1"\t"$3"\t"$2"\t"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$13"\t"$8}}' > windows.fa.DB.euX.bed


#################################################################################################
# Intersect euk and pro blast hits
#################################################################################################
# get overlapping HSPs between pro (tmp1) and euk (tmp2).
## for each pro hit, all overlapping euk hits are returned (bedtools intersect -a windows.fa.DB.pro.bed -b windows.fa.DB.euk.bed -wao)
### sort overlap by euk bitscore (sort -t$'\t' -k 1,4 -k 10,10gr)
### keep only the best (-u) scoring euk HSP (as sorted before) (-k10,10gr) for each pro HSP ( sort -t$'\t' -u -k1,4)
#################################################################################################

## blastn
bedtools intersect -a windows.fa.DB.pro.bed -b windows.fa.DB.euk.bed -wao |sort -t$'\t' -k 1,4 -k 10,10gr | sort -t$'\t' -u -k1,4  > windows.fa.DB.proVSeuk.bed

## blastx
bedtools intersect -a windows.fa.DB.prX.bed -b windows.fa.DB.euX.bed -wao |sort -t$'\t' -k 1,4 -k 10,10gr | sort -t$'\t' -u -k1,4  > windows.fa.DB.proVSeuk.X.bed


#################################################################################################
# Filter 1
#################################################################################################
# Filter hsps where bacterial bitscore is higher than eukaryotic and where bacterial bitscore is > 69.
## bitscore 69 ~ evalue 1e-5
## bitscore 73 ~ evalue 1e-10
#################################################################################################

##blastn
cat windows.fa.DB.proVSeuk.bed   |awk -F $'\t' '{if ($5>$10 && $5 > 69) print $0}'|sort -t$'\t' -k 1,3 -k 5,5gr | sort -t$'\t' -u -k1,3 > windows.fa.DB.LGTs.bed

##blastx
cat windows.fa.DB.proVSeuk.X.bed |awk -F $'\t' '{if ($5>$10 && $5 > 69) print $0}'|sort -t$'\t' -k 1,3 -k 5,5gr | sort -t$'\t' -u -k1,3 > windows.fa.DB.LGTs.X.bed

#################################################################################################
# Filter 2
#################################################################################################
# Filter LGT candidate windows
## keep only those where
### scaffold is not among the "contaminants"
### the bitscore difference is > 100
### the prokaryotic HSP is larger than 100 bp
#################################################################################################

## prepare grep-able list of contaminating scaffolds
cat contaminantscaffolds.csv|sed 1d|cut -f2|perl -pe 's/\n/\:\n/g' > contaminants

##blastn
### strict filter: pro HSP length > 100 bp  & bitscore diff > 100
bedtools sort -i windows.fa.DB.LGTs.bed  |sort -k 1,4 -u|grep -f contaminants - -v|awk -F $'\t' 'BEGIN {OFS = FS} {print $0"\t"$5-$10}'|sort -k 12,12gr|awk -F $'\t' 'BEGIN {OFS = FS} {if ($12 > 100 && sqrt(($2-$3)^2) > 100) print $0}' > windows.fa.DB.LGTs.filtered.bed
### loose filter: pro HSP length > 50 bp  & bitscore diff > 50
bedtools sort -i windows.fa.DB.LGTs.bed  |sort -k 1,4 -u|grep -f contaminants - -v|awk -F $'\t' 'BEGIN {OFS = FS} {print $0"\t"$5-$10}'|sort -k 12,12gr|awk -F $'\t' 'BEGIN {OFS = FS} {if ($12 > 50 && sqrt(($2-$3)^2) > 50) print $0}' > windows.fa.DB.LGTs.filtered.loose.bed
##blastx
bedtools sort -i windows.fa.DB.LGTs.X.bed|sort -k 1,4 -u|grep -f contaminants - -v|awk -F $'\t' 'BEGIN {OFS = FS} {print $0"\t"$5-$10}'|sort -k 12,12gr|awk -F $'\t' 'BEGIN {OFS = FS} {if ($12 > 100 && sqrt(($2-$3)^2) > 100) print $0}' > windows.fa.DB.LGTs.filtered.X.bed

## confirm contaminants
bedtools sort -i windows.fa.DB.LGTs.bed  |sort -k 1,4 -u|grep -f contaminants - |awk -F $'\t' 'BEGIN {OFS = FS} {print $0"\t"$5-$10}'|sort -k 12,12gr|awk -F $'\t' 'BEGIN {OFS = FS} {if ($12 < 0 ) print $0}' > windows.fa.DB.contaminants.filtered.bed
#NW_021768749.1 Nvit

#################################################################################################
# merge HSPs to loci
#################################################################################################
# split window coordinates in tabs (perl -pe 's/(^scaffold[0-9]+)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g')
# generate genomic position (awk  -F $'\t' 'BEGIN {OFS = FS} {print $1,$2+$4,$2+$5-1,$6,$7,$12,$14,$13}')
# merge loci <500 bp apart (bedtools merge -i - -c 4,5,6,7,8 -o first,collapse,collapse,collapse,collapse -d 500)
## report first pro hit (col 4)
## collapse bac bitscores (col 5)
## collapse euk bitscores (col 6)
## collapse pairwise differences in bitscores (col 7)
## collapse overlap (col 8)
#################################################################################################

# colnames of output
# scaffold  start end bestProHit;tstart;tend;evalue;bits;alnlen;pident;taxname; collapsed(proBit)  collapsed(eukBit)  collapsed(bitDiff) collapsed(overlapping)

##blastn
cat windows.fa.DB.LGTs.filtered.bed      | perl -pe 's/(^.+?)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g'|awk  -F $'\t' 'BEGIN {OFS = FS} {print $1,$2+$4,$2+$5-1,$6,$7,$12,$14,$13}'|bedtools sort -i - |bedtools merge -i - -c 4,5,6,7,8 -o first,collapse,collapse,collapse,collapse -d 500 |bedtools sort -i - > LGTs.candidateloci.bed
cat windows.fa.DB.LGTs.filtered.bed      | perl -pe 's/(^scaffold[0-9]+)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g'|awk  -F $'\t' 'BEGIN {OFS = FS} {print $1,$2+$4,$2+$5-1,$6,$7,$12,$14,$13}'|bedtools sort -i - |bedtools merge -i - -c 4,5,6,7,8 -o first,collapse,collapse,collapse,collapse -d 500 |bedtools sort -i - > LGTs.candidateloci.bed
cat windows.fa.DB.LGTs.filtered.loose.bed| perl -pe 's/(^scaffold[0-9]+)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g'|awk  -F $'\t' 'BEGIN {OFS = FS} {print $1,$2+$4,$2+$5-1,$6,$7,$12,$14,$13}'|bedtools sort -i - |bedtools merge -i - -c 4,5,6,7,8 -o first,collapse,collapse,collapse,collapse -d 100 |bedtools sort -i - > LGTs.candidateloci.loose.bed

## limit loose candidates to loci > 100 bp
cat LGTs.candidateloci.loose.bed |awk '{if (sqrt(($2-$3)^2) > 100 ) print $0}'

##blastx
cat windows.fa.DB.LGTs.filtered.X.bed    | perl -pe 's/(^scaffold[0-9]+)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g'|awk  -F $'\t' 'BEGIN {OFS = FS} {print $1,$2+$4,$2+$5-1,$6,$7,$12,$14,$13}'|bedtools sort -i - |bedtools merge -i - -c 4,5,6,7,8 -o first,collapse,collapse,collapse,collapse -d 500 > LGTs.candidateloci.X.bed

# Overlap blastx and blastn loci
#bedtools intersect -b LGTs.candidateloci.X.bed -a LGTs.candidateloci.bed -wao

#################################################################################################
# Analyses
#################################################################################################

## extract regions to fasta files
#cat ${query}_nextpolish.fasta|cut -f 1 -d "|" > genome.fa

bedtools getfasta -fi genome.fa -bed LGTs.candidateloci.loose.bed -fo LGTs.candidateloci.loose.fa
bedtools getfasta -fi genome.fa -bed LGTs.candidateloci.bed -fo LGTs.candidateloci.fa
bedtools getfasta -fi genome.fa -bed LGTs.candidateloci.X.bed -fo LGTs.candidateloci.X.fa

## retrieve overlapping blastx hits
cat windows.fa.DB.prX.bed | perl -pe 's/(^scaffold[0-9]+)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g'|awk  -F $'\t' 'BEGIN {OFS = FS} {print $1,$2+$4,$2+$5-1,$6}' > loci.DB.prX.bed
cat loci.DB.prX.bed| sort -t$'\t' -u -k1,3 --merge|bedtools sort -i -|bedtools merge -i - -c 4 -o collapse|bedtools intersect -a LGTs.candidateloci.bed -b - -wao|bedtools sort -i > LGTs.candidateloci.proteins.bed

## retrieve coverage in surrounding of loci
### plus/minus 8kb
cut -f 1,2 genome.fa.fai > genome.file
bedtools slop -b 8000 -i LGTs.candidateloci.bed -g genome.file|bedtools intersect -b genome.overlappingwindows.cov.tsv -a - -wao |bedtools merge -i - -c 4,12,13 -o first,collapse,collapse|bedtools sort -i > LGTs.candidateloci.coverage.bed
bedtools slop -b 8000 -i LGTs.candidateloci.X.bed -g genome.file|bedtools intersect -b genome.overlappingwindows.cov.tsv -a - -wao |bedtools merge -i - -c 4,12,13 -o first,collapse,collapse|bedtools sort -i > LGTs.candidateloci.X.coverage.bed


grep ntJoin19 genome.overlappingwindows.cov.tsv
bedtools intersect -b genome.overlappingwindows.cov.tsv -a LGTs.candidateloci.bed -wao |bedtools merge -i - -c 4,12,13 -o first,collapse,collapse|bedtools sort -i > LGTs.candidateloci.coverage.bed

# merge loci <5 kb apart
## blastn
bedtools merge -d 5000 -i LGTs.candidateloci.bed > LGTs.5kb.candidateregions.bed

#################################################################################################
# Unused
#################################################################################################

# ## implement screen for low complexity
# # dustmasker maybe not working, reports multiple lines per entry?
# dustmasker -in LGTs.candidateloci.fa -outfmt acclist > f1
# seqkit fx2tab LGTs.candidateloci.fa -g -G -l > f2
# awk 'NR==FNR {h[$1] = $0; next} {if (h[">"$1]!=""){print $0"\t"h[">"$1]}else{print $0"\t"$1"\t0\t0"}}' f1 f2 |awk '{print $1"\t"$3"\t"($8-$7)/$3}'|perl -pe 's/(^scaffold[0-9]+)\:([0-9]+)-([0-9]+)/$1\t$2\t$3/g'|bedtools sort -i - > LGTs.candidateloci.lowComplexity.bed
# bedtools intersect -b LGTs.candidateloci.lowComplexity.bed -a LGTs.candidateloci.coverage.bed -wao|cut -f 1-6,11 > LGTs.candidateloci.cov.lowComp.bed
#
# # preliminary blast analysis
# cat LGTs.candidateloci.bed |bedtools getfasta -bed - -fi genome.fa|blastx -query -  -db /global/databases/swissprot/2020_11/uniprot_sprot.fasta -outfmt 6 -num_threads 6 -evalue 1e-1 > LGTs.candidateloci.vsUniProt.bls
# cat LGTs.candidateloci.vsUniProt.bls|sort -k1,1 -k11,11g -k12,12gr | sort -u -k1,1 --merge > LGTs.candidateloci.vsUniProt.bh
#
# query=$(cut -f 2 -d "|" LGTs.candidateloci.vsUniProt.bh|uniq|tr "\n" "|"|perl -pe 's/\|$//g')
# export NCBI_API_KEY='7820ce7ebdfc0ed36d2d1d2fc47746f0c209'
# esearch -db protein -query $query|elink -target gene|efetch -format tabular

#################################################################################################
# Notes
#################################################################################################

## Known contaminations of the eukaryotic scaffold
## grep NW_021768749.1 windows.fa.DB.euk.bed
## gb|SSSL01000410.1|
## JACTOK010014482.1
## JACTOK010038917.1
## JXPW01070074.1
## JXPW01070074.1
## NW_021768749.1
## NW_023337194.1
## get all windows that hit good against bacteria and good against insects
cat windows.fa.DB.proVSeuk.bed   |awk -F $'\t' '{if ($5>200 && $10>200) print $0}'|cut -f 9|cut -f 8 -d ";"|sort|uniq -c
