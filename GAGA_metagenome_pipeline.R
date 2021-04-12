#!/usr/bin/env Rscript

## This metagenome pipeline is designed to classify GAGA genome assembly scaffolds into tentative taxonomic origins.
## The focus groups are Insecta and Bacteria. 
## The processes used for classification are
## 1. mmseqs blastn vs prok database (NCBI: 1908 complete prokaryotic genomes)
## 2. mmseqs blastn vs insect database (NCBI: 43 high quality complete genomes)
## 3. mmseqs blastx vs prok database (Uniprot90/50 database - Bacteria part)
## 4. mmseqs blastx vs insect database (Uniprot90/50 database - Insecta part)
## 5. barrnap rRNA prediction
## 6. calculation of sequencing coverage of scaffolds from raw reads
## 7. calculation of GC content of scaffolds

## Load packages
suppressPackageStartupMessages({
library(cowplot)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(GenomicRanges)
library(ggfortify)
library(ggpubr)
library(DescTools)
library(factoextra)
library(ade4)
})

## Prepare files
print("Loading files...")
# GAGA ID and result folder location (can be passed from the commandline)
args <- commandArgs(trailingOnly = T)
if (length(args) == 0) { 
  stop("Need to provide GAGA-ID")
} else {
  warning("No result directory is provided, using current directory for outputs...")
  id <- args[1]
  folder <- paste(getwd(),"/",sep="")
} 

# rRNA annotations 
proRNA.tmp<-read.csv(paste(folder,"genome.pro_blastn.rRNA.windows.bed",sep=""),sep='\t',header=F,comment.char = "#")
eukRNA.tmp<-read.csv(paste(folder,"genome.euk_blastn.rRNA.windows.bed",sep=""),sep='\t',header=F,comment.char = "#")
eukRNA.tmp<-eukRNA.tmp[,c(1:3,6,7,8,9,10,12)]
proRNA.tmp<-proRNA.tmp[,c(1:3,6,7,8,9,10,12)]
colnames(eukRNA.tmp)<-c("scf","start","end","eukRNA.type","eukRNA.start","eukRNA.end","eukRNA.eval","eukRNA.strand","eukRNA.att")
colnames(proRNA.tmp)<-c("scf","start","end","type","proRNA.start","proRNA.end","proRNA.eval","proRNA.strand","proRNA.att")
eukRNA.tmp$qseqid<-paste(eukRNA.tmp$scf,":",eukRNA.tmp$start,"-",eukRNA.tmp$end,sep="")
proRNA.tmp$qseqid<-paste(proRNA.tmp$scf,":",proRNA.tmp$start,"-",proRNA.tmp$end,sep="")
eukRNA<-eukRNA.tmp %>% 
  dplyr::group_by(qseqid) %>% 
  dplyr::arrange(eukRNA.eval) %>% 
  dplyr::slice(1)
proRNA<-proRNA.tmp %>% 
  dplyr::group_by(qseqid) %>% 
  dplyr::arrange(proRNA.eval) %>% 
  dplyr::slice(1)

# mmseqs blastn hits files
euk<-read.csv(paste(folder,"query.", id, ".DB.euk_blastn.bh",sep=""),sep='\t',header=F)
pro<-read.csv(paste(folder,"query.", id, ".DB.pro_blastn.bh",sep=""),sep='\t',header=F)
noAntInsect<-read.csv(paste(folder,"query.", id, ".DB.euk_noAnts_blastn.bh",sep=""),sep='\t',header=F)
human<-read.csv(paste(folder,"query.", id, ".DB.human_blastn.bh",sep=""),sep='\t',header=F)

# mmseqs blastx hits files
euk.x<-read.csv(paste(folder,"query.", id, ".DB.euk_blastx.bh",sep=""),sep='\t',header=F)
pro.x<-read.csv(paste(folder,"query.", id, ".DB.pro_blastx.bh",sep=""),sep='\t',header=F)

#  GC content & window Count per scaffold
GC<-read.csv(paste(folder,"genome.GC.tsv",sep=""),sep='\t',header=T)
wGC<-read.csv(paste(folder,"windows.GC.tsv",sep=""),sep='\t',header=T)
colnames(wGC)<-paste("w",colnames(wGC),sep="")
colnames(wGC)[1]<-"qseqid"
windowCount<-read.csv(paste(folder,"genome.overlappingwindows.tsv",sep=""),sep=':',header=F)

#coverage
cov<-read.csv(paste(folder,"genome.overlappingwindows.cov.tsv",sep=""),sep='\t',header=F)
cov$qseqid<-paste(cov$V1,":",cov$V2,"-",cov$V3,sep="")
colnames(cov)[c(1,4)]<-c("scaffold","cov")

# name columns
print("Assemblying into data frames...")

blsNcolnames<-c("qScf","qScfstart","qScfend","qseqid","qstart","qend","sseqid","sstart","send","evalue","bitscore","length","pident","sgi","sacc","stitle")
colnames(euk)<-paste("euk.",blsNcolnames,sep="")
colnames(euk.x)<-paste("euk.x.",blsNcolnames,sep="")
colnames(pro)<-paste("pro.",blsNcolnames,sep="")
colnames(pro.x)<-paste("pro.x.",blsNcolnames,sep="")
colnames(GC)<-c("scaffold","Length","GC")
colnames(noAntInsect)<-paste("noAntInsect.",blsNcolnames,sep="")
colnames(human)<-paste("human.",blsNcolnames,sep="")

# clean up prokaryotic taxonomy information
pro$pro.stitle<-gsub("Candidatus ","",pro$pro.stitle)
pro$pro.tax <- with(pro, ifelse(grepl(".*g_(.*?);.*", pro.sgi), gsub(".*g_(.*?);.*","\\1",pro.sgi), 
                         ifelse(grepl(".*f_(.*?);.*", pro.sgi), gsub(".*f_(.*?);.*","\\1",pro.sgi),
                         ifelse(grepl(".*o_(.*?);.*", pro.sgi), gsub(".*o_(.*?);.*","\\1",pro.sgi),
                         ifelse(grepl(".*c_(.*?);.*", pro.sgi), gsub(".*c_(.*?);.*","\\1",pro.sgi),
                         ifelse(grepl(".*p_(.*?);.*", pro.sgi), gsub(".*p_(.*?);.*","\\1",pro.sgi),
                                                                gsub(".*d_(.*?);.*","\\1",pro.sgi)))))))
pro.x$pro.x.stitle<-gsub("Candidatus ","",pro.x$pro.x.stitle)
pro.x$pro.x.tax <- with(pro.x, ifelse(grepl(".*g_(.*?);.*", pro.x.sgi), gsub(".*g_(.*?);.*","\\1",pro.x.sgi), 
                               ifelse(grepl(".*f_(.*?);.*", pro.x.sgi), gsub(".*f_(.*?);.*","\\1",pro.x.sgi),
                               ifelse(grepl(".*o_(.*?);.*", pro.x.sgi), gsub(".*o_(.*?);.*","\\1",pro.x.sgi),
                               ifelse(grepl(".*c_(.*?);.*", pro.x.sgi), gsub(".*c_(.*?);.*","\\1",pro.x.sgi),
                               ifelse(grepl(".*p_(.*?);.*", pro.x.sgi), gsub(".*p_(.*?);.*","\\1",pro.x.sgi),
                                                                        gsub(".*d_(.*?);.*","\\1",pro.x.sgi)))))))

euk$euk.tax <- with(euk, ifelse(grepl(".*g_(.*?);.*", euk.sgi), gsub(".*g_(.*?);.*","\\1",euk.sgi), 
                         ifelse(grepl(".*f_(.*?);.*", euk.sgi), gsub(".*f_(.*?);.*","\\1",euk.sgi),
                         ifelse(grepl(".*o_(.*?);.*", euk.sgi), gsub(".*o_(.*?);.*","\\1",euk.sgi),
                         ifelse(grepl(".*c_(.*?);.*", euk.sgi), gsub(".*c_(.*?);.*","\\1",euk.sgi),
                         ifelse(grepl(".*p_(.*?);.*", euk.sgi), gsub(".*p_(.*?);.*","\\1",euk.sgi),
                                                                gsub(".*d_(.*?);.*","\\1",euk.sgi)))))))

euk.x$euk.x.tax <- with(euk.x, ifelse(grepl(".*g_(.*?);.*", euk.x.sgi), gsub(".*g_(.*?);.*","\\1",euk.x.sgi), 
                               ifelse(grepl(".*f_(.*?);.*", euk.x.sgi), gsub(".*f_(.*?);.*","\\1",euk.x.sgi),
                               ifelse(grepl(".*o_(.*?);.*", euk.x.sgi), gsub(".*o_(.*?);.*","\\1",euk.x.sgi),
                               ifelse(grepl(".*c_(.*?);.*", euk.x.sgi), gsub(".*c_(.*?);.*","\\1",euk.x.sgi),
                               ifelse(grepl(".*p_(.*?);.*", euk.x.sgi), gsub(".*p_(.*?);.*","\\1",euk.x.sgi),
                                                                        gsub(".*d_(.*?);.*","\\1",euk.x.sgi)))))))

noAntInsect$noAntInsect.tax <- with(noAntInsect, ifelse(grepl(".*g_(.*?);.*", noAntInsect.sgi), gsub(".*g_(.*?);.*","\\1",noAntInsect.sgi), 
                                ifelse(grepl(".*f_(.*?);.*", noAntInsect.sgi), gsub(".*f_(.*?);.*","\\1",noAntInsect.sgi),
                                       ifelse(grepl(".*o_(.*?);.*", noAntInsect.sgi), gsub(".*o_(.*?);.*","\\1",noAntInsect.sgi),
                                              ifelse(grepl(".*c_(.*?);.*", noAntInsect.sgi), gsub(".*c_(.*?);.*","\\1",noAntInsect.sgi),
                                                     ifelse(grepl(".*p_(.*?);.*", noAntInsect.sgi), gsub(".*p_(.*?);.*","\\1",noAntInsect.sgi),
                                                            gsub(".*d_(.*?);.*","\\1",noAntInsect.sgi)))))))

human$human.tax <- with(human, ifelse(grepl(".*g_(.*?);.*", human.sgi), gsub(".*g_(.*?);.*","\\1",human.sgi), 
                                                        ifelse(grepl(".*f_(.*?);.*", human.sgi), gsub(".*f_(.*?);.*","\\1",human.sgi),
                                                               ifelse(grepl(".*o_(.*?);.*", human.sgi), gsub(".*o_(.*?);.*","\\1",human.sgi),
                                                                      ifelse(grepl(".*c_(.*?);.*", human.sgi), gsub(".*c_(.*?);.*","\\1",human.sgi),
                                                                             ifelse(grepl(".*p_(.*?);.*", human.sgi), gsub(".*p_(.*?);.*","\\1",human.sgi),
                                                                                    gsub(".*d_(.*?);.*","\\1",human.sgi)))))))

# merge pro and euk blastn/blastx results
m1<- merge(euk,pro,by.x="euk.qseqid",by.y="pro.qseqid",all.x=T,all.y=T) %>%
  merge(.,euk.x,by.x="euk.qseqid",by.y="euk.x.qseqid",all.x=T,all.y=T) %>%
  merge(.,pro.x,by.x="euk.qseqid",by.y="pro.x.qseqid",all.x=T,all.y=T) %>%
  merge(.,noAntInsect,by.x="euk.qseqid",by.y="noAntInsect.qseqid",all.x=T,all.y=T) %>%
  merge(.,human,by.x="euk.qseqid",by.y="human.qseqid",all.x=T,all.y=T) %>%
  merge(.,eukRNA,by.x="euk.qseqid",by.y="qseqid",all.x=T,all.y=T) %>%
  merge(.,proRNA,by.x="euk.qseqid",by.y="qseqid",all.x=T,all.y=T)

colnames(m1)[colnames(m1)=="euk.qseqid"]<-"qseqid"
m1$scaffold<-gsub(":.*","",m1$qseqid,perl=T)
m1$start<-as.numeric(gsub(".*:(.*)-.*","\\1",m1$qseqid,perl=T))
m1$end<-as.numeric(gsub(".*:.*-(.*)","\\1",m1$qseqid,perl=T))

#select relevant columns to keep
m2<-m1[,c("scaffold","start","end","qseqid",
          "euk.evalue","euk.bitscore","euk.pident","euk.sstart","euk.send","euk.qstart","euk.qend","euk.stitle","euk.tax",
          "pro.evalue","pro.bitscore","pro.pident","pro.sstart","pro.send","pro.qstart","pro.qend","pro.stitle","pro.tax",
          "noAntInsect.evalue","noAntInsect.bitscore","noAntInsect.pident","noAntInsect.sstart","noAntInsect.send","noAntInsect.stitle","noAntInsect.tax",
          "human.evalue","human.bitscore","human.pident","human.sstart","human.send","human.stitle","human.tax",
          "eukRNA.eval","proRNA.eval","euk.x.bitscore","euk.x.tax","pro.x.bitscore","pro.x.tax","pro.sacc")]
m2$scaffold<-as.factor(m2$scaffold)

#calculate coverage across scaffolds
print("Calculating coverage of scaffolds...")

cov$relCov<-log((cov$cov+1)/(median(cov$cov)),2)
covScf<-cov %>% 
  group_by(scaffold) %>%
  summarise(no_rows = mean(relCov))
colnames(covScf)<-c("scaffold","coverage")

# merge blastn/x results with GC content, coverage and windows-wise GC content
m4<-merge(m2,cov[,c("qseqid","relCov")],by="qseqid",all.x=T,all.y=T) %>% 
  merge(.,wGC[,c("qseqid","wX.GC")],by="qseqid",all.x=T) %>%
  merge(.,GC,by.x="scaffold",by.y="scaffold",all.x=T)
m4$scaffold<-gsub(":.*","",m4$qseqid,perl=T)

#order scaffolds
m4<-m4[order(m4$scaffold,m4$start,decreasing = F),]

# tag those windows that are rRNAs
m4$pro.tax[m4$eukRNA.eval>m4$proRNA.eval | (!is.na(m4$proRNA.eval) & is.na(m4$eukRNA.eval))]<-"pro.rRNA"
m4$pro.tax[m4$eukRNA.eval<m4$proRNA.eval | (is.na(m4$proRNA.eval) & !is.na(m4$eukRNA.eval))]<-"euk.rRNA"

# calculate pro-euk bitscore difference for each window
m4$window.bsdiff <- with(m4, ifelse(is.na(pro.bitscore) & !is.na(euk.bitscore), 0 - euk.bitscore,
                                    ifelse(is.na(euk.bitscore) & !is.na(pro.bitscore), pro.bitscore, 
                                           pro.bitscore - euk.bitscore)))

# calculate the longest consecutive window-best-hit for euk & pro along each entire scaffold 
longestConsecutiveWindows <- function(v) {
  v <- v[!is.na(v)]
  v <- ifelse(v>0, T, F)
  return(c(length(v[v==T]), length(v[v==F])))
}

scaffoldHitStats<-m4 %>% 
  dplyr::group_by(scaffold) %>%
  dplyr::summarise(windowCount = length(scaffold),
                   euk_windows = round(length(window.bsdiff[window.bsdiff < 0 & !is.na(window.bsdiff)])/length(scaffold),2),
                   pro_windows = round(length(window.bsdiff[window.bsdiff > 0 & !is.na(window.bsdiff)])/length(scaffold),2),
                   amb_windows = round(length(window.bsdiff[window.bsdiff == 0 & !is.na(window.bsdiff)])/length(scaffold),2),
                   na_windows = round(length(window.bsdiff[is.na(window.bsdiff)])/length(scaffold),2),
                   longestProWindowPerc = round(longestConsecutiveWindows(window.bsdiff)[1]/length(scaffold),2), 
                   longestProWindowSize = longestConsecutiveWindows(window.bsdiff)[1],
                   longestProWindowMeanCov = mean(relCov))
colnames(scaffoldHitStats)<-c("scaffold","window.count","euk.wind.perc","pro.wind.perc","amb.wind.perc",
                              "no-hit.wind.perc","LongestContProWindow.perc","LongestContProWindowSize",
                              "longestProWindowMeanCov")

######################## experimental steps #################################
# # calculate the number of switches of window-best-hit (euk<->pro) along each entire scaffold 
# # a function to count pro<->euk hit switch 
# hitSwitchCount <- function(v) {
#   v <- v[!is.na(v)]
#   v <- ifelse(v>0, T, F)
#   sum(diff(v) == 1 | diff(v) == -1)
# }
# # a function to calculate a score to capture hit patterns
# # if score is negative, the value closer to -1 represents smaller amount of pro-window and larger amount of euk-window
# # if score is negative, the value closer to 0 represents both smaller amount of pro-window and euk-window, lots of no-hit
# # if score is positive, the value closer to 1 represents larger amount of pro-window and smaller amount of euk-window
# # if score is positive, the value closer to 0 represents both smaller amount of pro-window and euk-window, lots of no-hit
# hitPatternScore <- function(switchCount, windowBSdiff, scaffoldLen) {
#   ProEukwindowDiff <- length(windowBSdiff[windowBSdiff > 0]) - length(windowBSdiff[windowBSdiff < 0])
#   switchCount * (ProEukwindowDiff / scaffoldLen)
# }
# windowHitStats<-m4 %>% 
#   dplyr::group_by(scaffold) %>%
#   dplyr::summarise(no_rows = length(scaffold), 
#                    euk_bs = sum(euk.bitscore,na.rm=T), 
#                    pro_bs = sum(pro.bitscore,na.rm=T), 
#                    euk_windows = length(window.bsdiff[window.bsdiff < 0]) / length(scaffold),
#                    pro_windows = length(window.bsdiff[window.bsdiff > 0]) / length(scaffold),
#                    ProEukSwitchCount = hitSwitchCount(window.bsdiff), 
#                    patternScore = hitPatternScore(hitSwitchCount(window.bsdiff), 
#                                                   window.bsdiff, length(scaffold)))
# colnames(windowHitStats)<-c("scaffold","WindowCount","euk.bitscore","pro.bitscore",
#                             "euk.hits.perc","pro.hits.perc","ProEukSwitchCount",
#                             "WindowSwitchPatternScore")
###############################################################################

# (blastn) Filter all windows with better hit against euk
print("Classifying sliding windows...")
eukWindows<-subset(m4,
                   (euk.bitscore>pro.bitscore  |   # is euk.bitscore higher than pro.bitscore
                   (eukRNA.eval<proRNA.eval & (is.na(euk.bitscore) & is.na(pro.bitscore)))  |   # is rRNA evalue higher agains pro than in euk
                   (!is.na(euk.bitscore) & is.na(pro.bitscore)) | # do we have no prok.bitscore but a euk.bitscore?
                   (!is.na(eukRNA.eval) & is.na(proRNA.eval) & (is.na(euk.bitscore) & is.na(pro.bitscore)))))    # do we have no prok.rRNA hit but a euk.rRNA hit?

# (blastx) Filter all windows with better hit against euk
eukWindows.x<-subset(m4,
                     euk.x.bitscore>pro.x.bitscore |   # is euk.x.bitscore higher than pro.x.bitscore
                    (!is.na(euk.x.bitscore) & is.na(pro.x.bitscore))) # do we have no prok.x.bitscore but a euk.x.bitscore?

# (blastn) Filter all windows with better hit against pro
proWindows<-subset(m4,   
                   (euk.bitscore<=pro.bitscore |   # is euk.bitscore lower than pro.bitscore
                   (eukRNA.eval>=proRNA.eval & (is.na(euk.bitscore) & is.na(pro.bitscore))) |   # is rRNA evalue lower agains pro than in euk
                   (is.na(euk.bitscore) & !is.na(pro.bitscore)) | # do we have a prok.bitscore but no euk.bitscore?
                   (is.na(eukRNA.eval) & !is.na(proRNA.eval) & (is.na(euk.bitscore) & is.na(pro.bitscore)))))  # do we have a prok.rRNA hit but no euk.rRNA hit?

# (blastx) Filter all windows with better hit against pro
proWindows.x<-subset(m4,
                     euk.x.bitscore<=pro.x.bitscore |   # is euk.x.bitscore lower than pro.x.bitscore
                    (is.na(euk.x.bitscore) & !is.na(pro.x.bitscore))) # do we have a prok.bitscore but no euk.bitscore?

# (blastn) Filter all windows with better hit against human
humanWindows<-subset(m4,   
                    (!is.na(human.bitscore) & is.na(pro.bitscore) & is.na(euk.bitscore)) | # do we only have a human.bitscore
                    (!is.na(human.bitscore) & !is.na(pro.bitscore) & !is.na(euk.bitscore) & 
                    human.bitscore > pro.bitscore & human.bitscore > euk.bitscore) |
                    (!is.na(human.bitscore) & !is.na(pro.bitscore) & is.na(euk.bitscore) &
                    human.bitscore > pro.bitscore) |
                    (!is.na(human.bitscore) & is.na(pro.bitscore) & !is.na(euk.bitscore) &
                    human.bitscore > euk.bitscore))

## Identify rate of "prokaryotic windows" across entire chromosom
print("Gathering sliding-window inforation into scafolds...")

## calculate scaffolds-wide information
# get number of windows producing significant hits vs insectDB
eukWindowCount<-eukWindows %>% 
  dplyr::group_by(scaffold) %>%
  dplyr::summarise(no_rows = length(scaffold),bs=sum(euk.bitscore,na.rm=T))
colnames(eukWindowCount)<-c("scaffold","eukWindows","euk.bitscore")

eukWindowCount.x<-eukWindows.x %>% 
  dplyr::group_by(scaffold) %>%
  dplyr::summarise(no_rows = length(scaffold),bs=sum(euk.x.bitscore,na.rm=T))
colnames(eukWindowCount.x)<-c("scaffold","eukWindows.x","euk.x.bitscore")

# get number of windows producing significant hits vs prokDB
proWindowCount<-proWindows %>% 
  dplyr::group_by(scaffold) %>%
  dplyr::summarise(no_rows = length(scaffold),bs=sum(pro.bitscore,na.rm=T))
colnames(proWindowCount)<-c("scaffold","proWindows","pro.bitscore")

proWindowCount.x<-proWindows.x %>% 
  dplyr::group_by(scaffold) %>%
  dplyr::summarise(no_rows = length(scaffold),bs=sum(pro.x.bitscore,na.rm=T))
colnames(proWindowCount.x)<-c("scaffold","proWindows.x","pro.x.bitscore")

humanWindowCount<-humanWindows %>% 
  dplyr::group_by(scaffold) %>%
  dplyr::summarise(no_rows = length(scaffold),bs=sum(human.bitscore,na.rm=T))
colnames(humanWindowCount)<-c("scaffold","humanWindows","human.bitscore")

chrSummary<-merge(GC,proWindowCount,by.x="scaffold",by.y="scaffold",all.y=T,all.x=T) %>% 
  merge(.,proWindowCount.x,by.x="scaffold",by.y="scaffold",all.y=T,all.x=T) %>%
  merge(.,eukWindowCount,by.x="scaffold",by.y="scaffold",all.x=T,all.y=T) %>%
  merge(.,eukWindowCount.x,by.x="scaffold",by.y="scaffold",all.x=T,all.y=T) %>%
  merge(.,humanWindowCount,by.x="scaffold",by.y="scaffold",all.x=T,all.y=T)

# set NAs to 0
chrSummary[is.na(chrSummary)]<-0 
# calculate ratio of windows that have a better hit against the prokaryotic DB than against the eukaryotic DB
chrSummary$ratio<-chrSummary$proWindows/(chrSummary$proWindows+chrSummary$eukWindows)
chrSummary$bsratio<-chrSummary$pro.bitscore-chrSummary$euk.bitscore
chrSummary$ratio.x<-chrSummary$proWindows.x/(chrSummary$proWindows.x+chrSummary$eukWindows.x)
chrSummary$bsratio.x<-chrSummary$pro.x.bitscore-chrSummary$euk.x.bitscore
chrSummary[is.na(chrSummary)]<-0

### Retrieve most frequently hit pro_taxon for each scaffold
proWindowsList<-split(proWindows,f = proWindows$scaffold,drop = T)

pro.tmp<-lapply(proWindowsList, FUN= function(x) {  
  x %>% 
    dplyr::group_by(pro.tax) %>%
    dplyr::summarise(count = length(pro.tax)) %>%
    dplyr::arrange(desc(count)) %>%
    dplyr::slice(1)
})

pro.bestMatch<-plyr::ldply(pro.tmp, rbind)

## Retrieve most frequently hit euk_taxon for each scaffold
eukWindowsList<-split(eukWindows,f = eukWindows$scaffold,drop = T)

euk.tmp<-lapply(eukWindowsList, FUN= function(x) {  
  x %>% 
    dplyr::group_by(euk.tax) %>%
    dplyr::summarise(count = length(euk.tax)) %>%
    dplyr::arrange(desc(count)) %>%
    dplyr::slice(1)
})

euk.bestMatch<-plyr::ldply(euk.tmp, rbind)

# Stop the analysis if there's zero pro-window detected
if (nrow(proWindows) == 0) { stop("No bacterial window detected! Analysis stops.") }

#Combine all scaffold-wide information 
chrSum<-merge(chrSummary,euk.bestMatch[,1:2],by.x="scaffold",by.y=".id",all.x=T) %>%
  merge(.,pro.bestMatch[,1:2],by.x="scaffold",by.y=".id",all.x=T) %>%
  merge(.,covScf,by="scaffold",all.x=T) %>%
  merge(.,scaffoldHitStats,by="scaffold",all.x=T)

# classify a scaffold as prokaryotic if 
#   1. 100 % of hitted windows are identified as prokaryotic OR
#   2. more than 50% of hitted windows are identified as prokaryotic AND bsratio > 200
#   3. more than 50% of hitted windows are identified as prokaryotic AND 0 < bsratio < 200 AND thier
#     3.1 bsratio.x > 50 AND
#     3.2 GC content does not lie within 95% confidence interval of euk scaffolds GC content distribution AND
#     3.3 coverage does not lie within 95% confidence interval of euk scaffolds coverage distribution 
# note that it is possible to have ratio 1 for only small amount of proWindow but the rest are uncertain (no hit)
# for downstream genome annotation, recommending only the scaffolds tagged "pro" should be excluded
# for downstram microbiome analyses, recommending to include all "pro" and "uncertain" scaffolds
# Checking for human contamination has added: kingdom will be assigned human if more than half of
# the windows hitted better against human DB

print("Classifying scaffolds...")

chrSum.euk <- chrSum %>% filter(ratio == 0)
GC_medianCI.euk <- MedianCI(chrSum.euk$GC, method = "boot", conf.level = 0.95)
Coverage_medianCI.euk <- MedianCI(chrSum.euk$coverage, method = "boot", conf.level = 0.95)

chrSum$kingdom<-ifelse(chrSum$humanWindows>=chrSum$window.count/2,"human",
                ifelse(chrSum$ratio==1, "pro",
                ifelse(chrSum$ratio<=0.5, "euk",
                ifelse(chrSum$bsratio>200, "pro",
                ifelse(chrSum$bsratio>0 & chrSum$bsratio.x>=50 & chrSum$GC<GC_medianCI.euk[2] & chrSum$GC>GC_medianCI.euk[3] &
                       chrSum$coverage<Coverage_medianCI.euk[2] & chrSum$coverage>Coverage_medianCI.euk[3],"pro",
                ifelse(is.na(chrSum$ratio) & is.na(chrSum$bsratio.x),"uncertain","unknown"))))))

chrSum$type<-ifelse(chrSum$kingdom=="unknown", "unknown", 
             ifelse(chrSum$kingdom=="uncertain", "uncertain", 
             ifelse(chrSum$kingdom=="human", "human",
             ifelse(chrSum$kingdom=="pro", chrSum$pro.tax, 
             ifelse(chrSum$kingdom=="euk", chrSum$euk.tax, "ERROR")))))

chrSum$type[is.na(chrSum$type)] <- "unknown"

chrSum$kingdom.bitscore<-ifelse(chrSum$kingdom=="pro", chrSum$pro.bitscore,
                         ifelse(chrSum$kingdom=="human", chrSum$human.bitscore,
                         ifelse(chrSum$kingdom=="euk", chrSum$euk.bitscore, 0)))

# flag scaffolds where the cumulative bitscore is below 200 ("weak evidence")
chrSum$Evidence<-ifelse(chrSum$kingdom.bitscore<200 | chrSum$type=="unknown" | chrSum$type=="uncertain", "weak", "strong")

## generate an overview of contaminant scaffolds
print("Generating overviews of contaminant scaffolds...")
# human contaminants
contaminants.human.relax <- subset(chrSum,humanWindows>0)
contaminants.human <- subset(chrSum,kingdom=="human")

# a more relaxed criterion: all scaffolds with any bacterial window
contaminants.relax <- subset(chrSum,proWindows>0)
contaminantScaffoldSummary.relax<-contaminants.relax[, c("scaffold","kingdom","pro.tax","proWindows","eukWindows",
                                                         "pro.wind.perc","euk.wind.perc","amb.wind.perc",
                                                         "no-hit.wind.perc","LongestContProWindowSize",
                                                         "longestProWindowMeanCov","kingdom.bitscore",
                                                         "ratio","bsratio","ratio.x","bsratio.x",
                                                         "Evidence","window.count","Length","GC")]
contaminantScaffoldSummary.relax<-contaminantScaffoldSummary.relax[order(contaminantScaffoldSummary.relax$proWindows,
                                                                         decreasing = T),]

# a more stringent criterion: all scaffolds that meet the original selection criteria above
contaminants<-subset(chrSum,kingdom=="pro")

# generate contaminant scaffold summary table
contaminantScaffoldSummary<-contaminants[, c("scaffold","pro.tax","proWindows","eukWindows","pro.wind.perc",
                                             "euk.wind.perc","amb.wind.perc","no-hit.wind.perc",
                                             "LongestContProWindowSize","longestProWindowMeanCov",
                                             "kingdom.bitscore","ratio","bsratio","ratio.x","bsratio.x",
                                             "Evidence","window.count","Length","GC")]
contaminantScaffoldSummary<-contaminantScaffoldSummary[order(contaminantScaffoldSummary$proWindows,decreasing = T),]

# generate contaminant windows detail table
contaminantWindows.relax<-merge(proWindows[proWindows$scaffold %in% contaminantScaffoldSummary.relax$scaffold,],
                                eukWindows[eukWindows$scaffold %in% contaminantScaffoldSummary.relax$scaffold,],
                                by="qseqid",all.x=T)
contaminantWindows<-proWindows[proWindows$scaffold %in% contaminantScaffoldSummary$scaffold,]

# generate a simple contaminant table (each raw represents a bac taxon)
contaminantsTable.relax<-contaminants.relax %>% 
  group_by(pro.tax) %>%
  summarise(coverage = round(mean(coverage),3),gc=round(mean(GC),3),scaffolds=n())
contaminantsTable.relax$prevalence <- with(contaminantsTable.relax, scaffolds/nrow(chrSum))
contaminantsTable.relax<-contaminantsTable.relax[order(contaminantsTable.relax$scaffolds,decreasing = T),]

contaminantsTable<-contaminants %>% 
                   group_by(pro.tax) %>%
                   summarise(coverage = round(mean(coverage),3),gc=round(mean(GC),3),scaffolds=n())
contaminantsTable$prevalence <- with(contaminantsTable, scaffolds/nrow(chrSum))
contaminantsTable<-contaminantsTable[order(contaminantsTable$scaffolds,decreasing = T),]

OverviewTable <- matrix(c(nrow(contaminants),sum(contaminants$Length),
                          sum(contaminants$Length)/sum(chrSum$Length),
                          contaminantsTable$pro.tax[1],
                          length(contaminantScaffoldSummary.relax[contaminantScaffoldSummary.relax$LongestContProWindowSize>=20,]$scaffold),
                          length(contaminantScaffoldSummary.relax[contaminantScaffoldSummary.relax$LongestContProWindowSize>=20 & contaminantScaffoldSummary.relax$Length>=100000,]$scaffold),
                          nrow(contaminants.human)),ncol=7)
colnames(OverviewTable) <- c("Pro scaffolds count","Pro scaffolds length",
                             "Pro scaffolds total prevalance",
                             "Top pro tax","Putative misassemblies count",
                             "Putative misassemblies (>100kb) count",
                             "human scaffold count")

# save/print Table for the identified contaminants
print("Generating bacterial scaffolds list/report...")
write.table(OverviewTable,file=paste(folder,"OverviewTable.csv",sep=""),quote=F,sep ="\t",row.names=F)

write.table(contaminantsTable,file=paste(folder,"contaminantsTable.csv",sep=""),quote=F,sep="\t",row.names=F)
write.table(contaminantScaffoldSummary,file=paste(folder,"contaminantscaffolds.csv",sep=""),quote=F,sep="\t",row.names=F)
write.table(contaminantWindows,file=paste(folder,"contaminantwindows.csv",sep=""),quote=F,sep ="\t",row.names=F)

write.table(contaminantsTable.relax,file=paste(folder,"contaminantsTable.relax.csv",sep=""),quote=F,sep="\t",row.names=F)
write.table(contaminantScaffoldSummary.relax,file=paste(folder,"contaminantscaffolds.relax.csv",sep=""),quote=F,sep="\t",row.names=F)
write.table(contaminantWindows.relax,file=paste(folder,"contaminantwindows.relax.csv",sep=""),quote=F,sep="\t",row.names=F)

write.table(contaminants.human,file=paste(folder,"contaminants.human.csv",sep=""),quote=F,sep="\t",row.names=F)
write.table(contaminants.human.relax,file=paste(folder,"contaminants.human.relax.csv",sep=""),quote=F,sep="\t",row.names=F)

## extract top XX euk/pro taxa (to simplify the plot legend)
type.euk.df <- as.data.frame(table(as.factor(chrSum$type[chrSum$kingdom=="euk"])))
type.euk.df.ordered <- type.euk.df[order(-type.euk.df$Freq),]
top10_type.euk.df <- type.euk.df.ordered[1:10,]
top10_type.euk.char <- as.character(top10_type.euk.df$Var1)
top10_type.euk.char <- append(top10_type.euk.char, "other sp.")
top1_type.euk.char <- as.character(top10_type.euk.df$Var1)[1]
top1_type.euk.char <- append(top1_type.euk.char, "other euk sp.")

type.pro.df <- as.data.frame(table(as.factor(chrSum$type[chrSum$kingdom=="pro"])))
if (length(type.pro.df[,]) > 0) {
  type.pro.df.ordered <- type.pro.df[order(-type.pro.df$Freq),]
  top_type.pro.df <- ifelse(length(type.pro.df.ordered[,1]) >= 10, 
                              type.pro.df.ordered[1:10,], type.pro.df.ordered[,])
  top_type.pro.char <- as.character(top_type.pro.df$Var1)
  top_type.pro.char <- append(top_type.pro.char, "other pro sp.")
  top_fewer_type.pro.char <- ifelse(length(top_type.pro.df$Var1) >= 5,
                               as.character(top_type.pro.df$Var1)[1:5],
                               as.character(top_type.pro.df$Var1))
  top_fewer_type.pro.char <- append(top_fewer_type.pro.char, "other pro sp.")

  chrSum$plot_type <- ifelse(chrSum$kingdom=="pro",
                             top_fewer_type.pro.char[match(chrSum$type, top_fewer_type.pro.char, 
                                                  nomatch = tail(top_fewer_type.pro.char,n=1))],
                             ifelse(chrSum$kingdom=="euk",
                                  top1_type.euk.char[match(chrSum$type, top1_type.euk.char,
                                                  nomatch = 2)],chrSum$kingdom))
# order the legend display to show top 1 euk first then top pro in decrease freq.
  chrSum$plot_type <- factor(chrSum$plot_type, levels=c(top1_type.euk.char,top_fewer_type.pro.char))
}

## Plot pro window percentage vs GC
# kingdom (euk-or-pro) plot
print("Generating some plot...")
pmain<-ggplot(chrSum, aes(x=GC, y=coverage)) + 
  geom_point(aes(size=Length/1e+6,fill=kingdom,color=Evidence),pch=21,
             alpha=ifelse(chrSum$Evidence=="weak",0.1,.7),stroke=.7)+
  theme_classic()+ggtitle(id)+
  ylab(label = "relative Coverage (log2)")+
  scale_color_manual(name="Evidence",values=c("black","white"))+
  scale_size(name="Size (Mbp)")+
  scale_fill_discrete(name="")+
  theme(legend.spacing = unit(.05,"mm"))
ggsave(paste(folder,"Taxa_screen.kingdom.pdf",sep=""), pmain, device = "pdf")

# top pro taxa plot
if (length(type.pro.df[,]) > 0) {
  pmain.top5<-ggplot(chrSum, aes(x=GC, y=coverage)) + 
  geom_point(aes(size=Length/1e+6,fill=plot_type,color=Evidence),pch=21,
             alpha=ifelse(chrSum$Evidence=="weak",0.1,.7),stroke=.7)+
  theme_classic()+ggtitle(id)+
  ylab(label = "relative Coverage (log2)")+
  scale_color_manual(name="Evidence",values=c("black","white"))+
  scale_size(name="Size (Mbp)")+
  scale_fill_brewer(palette="Set1",name="Top freq. taxa")+
  theme(legend.spacing = unit(.05,"mm"))
  ggsave(paste(folder,"Taxa_screen.top5protaxa.pdf",sep=""), pmain.top5, device = "pdf")
}
# PCA plot
chrSum.pca_df <- chrSum[c("scaffold", "GC", "ratio", "coverage", "kingdom", "Evidence")] %>% 
  tibble::column_to_rownames(var = "scaffold") %>%
  na.omit()

chrSum.pca_df.pca <- dudi.pca(chrSum.pca_df[1:3], nf=3, scannf = FALSE)
p.pca <- fviz_pca_biplot(chrSum.pca_df.pca, label = c("var", "quali"), 
                         fill.ind = chrSum.pca_df$kingdom, col.ind = chrSum.pca_df$Evidence, 
                         pointshape = 21, palette = c("red","blue","white"),
                         title = id, legend.title = list(color = "Evidence", fill = "Taxa"), mean.point = F 
)

ggsave(paste(folder,"Taxa_screen.pca.pdf",sep=""), p.pca, device = "pdf")

# Barplot of contaminant scaffolds pro-/euk-window distribution
contaminants.hist.df <- contaminantScaffoldSummary.relax %>% 
                        select(scaffold,kingdom,proWindows,eukWindows,window.count,
                               LongestContProWindowSize,longestProWindowMeanCov) %>% 
                        mutate(
                          scaffold = gsub(".*scaffold[^[:digit:]]*(\\d+)","\\1",scaffold,ignore.case=T),
                          proWindows = ifelse(proWindows>0, log(proWindows), 0),
                          eukWindows = ifelse(eukWindows>0, -log(eukWindows), 0),
                          LongestContProWindowSize = ifelse(LongestContProWindowSize>0, 
                                                            log(LongestContProWindowSize), 0)
                        )
hist.p1 <- ggplot(contaminants.hist.df,aes(x=scaffold)) + 
  geom_bar(aes(y=proWindows,fill=kingdom),stat="identity") + 
  geom_bar(aes(y=eukWindows,fill=kingdom),stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  geom_hline(aes(yintercept=0), color="black") + 
  geom_point(aes(y=LongestContProWindowSize,shape="longest wind."),
             stat="identity",color="blue",size=2,alpha=.7) + 
  geom_line(aes(y=longestProWindowMeanCov,group=1,color="cyan3")) + 
  scale_color_identity(name=NULL,labels="rel. cov.",guide="legend") + 
  annotate(geom="text",x=-Inf,y=-Inf,label="eukWindow",hjust=0,vjust=-1) + 
  annotate(geom="text",x=-Inf,y=Inf,label="proWindow",hjust=0,vjust=1) +
  labs(y="Window count",shape=NULL) + ggtitle(id)

ggsave(paste(folder,"Scaffold_hits.hist.pdf",sep=""), hist.p1, device = "pdf")

