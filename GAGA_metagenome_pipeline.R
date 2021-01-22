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

# mmseqs blastx hits files
euk.x<-read.csv(paste(folder,"query.", id, ".DB.euk_blastx.bh",sep=""),sep='\t',header=F)
pro.x<-read.csv(paste(folder,"query.", id, ".DB.pro_blastx.bh",sep=""),sep='\t',header=F)

#  GC content & window Count per scaffold
GC<-read.csv(paste(folder,"genome.GC.tsv",sep=""),sep='\t',header=T)
wGC<-read.csv(paste(folder,"windows.GC.tsv",sep=""),sep='\t',header=T)
colnames(wGC)<-paste("w",colnames(wGC),sep="")
colnames(wGC)[1]<-"qseqid"
windowCount<-read.csv(paste(folder,"genome.overlappingwindows.tsv",sep=""),sep='\t',header=F)

#coverage
cov<-read.csv(paste(folder,"genome.overlappingwindows.cov.tsv",sep=""),sep='\t',header=F)
cov$qseqid<-paste(cov$V1,":",cov$V2,"-",cov$V3,sep="")
colnames(cov)[4]<-"cov"

# name columns
print("Assemblying into data frames...")

blsNcolnames<-c("qScf","qScfstart","qScfend","qseqid","qstart","qend","sseqid","sstart","send","evalue","bitscore","length","pident","sgi","sacc","stitle")
colnames(euk)<-paste("euk.",blsNcolnames,sep="")
colnames(euk.x)<-paste("euk.x.",blsNcolnames,sep="")
colnames(pro)<-paste("pro.",blsNcolnames,sep="")
colnames(pro.x)<-paste("pro.x.",blsNcolnames,sep="")
colnames(GC)<-c("scaffold","Length","GC")

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

# merge pro and euk blastn/blastx results
m1<- merge(euk,pro,by.x="euk.qseqid",by.y="pro.qseqid",all.x=T,all.y=T) %>%
  merge(.,euk.x,by.x="euk.qseqid",by.y="euk.x.qseqid",all.x=T,all.y=T) %>%
  merge(.,pro.x,by.x="euk.qseqid",by.y="pro.x.qseqid",all.x=T,all.y=T) %>%
  merge(.,eukRNA,by.x="euk.qseqid",by.y="qseqid",all.x=T,all.y=T) %>%
  merge(.,proRNA,by.x="euk.qseqid",by.y="qseqid",all.x=T,all.y=T)

colnames(m1)[colnames(m1)=="euk.qseqid"]<-"qseqid"
m1$scaffold<-gsub(":.*","",m1$qseqid,perl=T)
m1$start<-as.numeric(gsub(".*:(.*)-.*","\\1",m1$qseqid,perl=T))
m1$end<-as.numeric(gsub(".*:.*-(.*)","\\1",m1$qseqid,perl=T))

#select relevant columns to keep
m2<-m1[,c("scaffold","start","end", "qseqid","euk.evalue","euk.bitscore","euk.pident",
          "euk.sstart","euk.send","pro.evalue","pro.bitscore","pro.pident","euk.stitle",
          "pro.stitle","pro.sstart","pro.send","pro.tax","euk.tax","eukRNA.eval",
          "euk.x.bitscore","pro.x.bitscore","euk.x.tax","pro.x.tax",
          "proRNA.eval","euk.qstart","euk.qend","pro.qstart","pro.qend","pro.sacc")]
m2$scaffold<-as.factor(m2$scaffold)

# merge blastn/x results with GC content, coverage and windows-wise GC content
m4<-merge(m2,GC,by.x="scaffold",by.y="scaffold",all.x=T,all.y=T) %>% 
  merge(.,cov[,c("qseqid","cov")],by="qseqid",all.x=T) %>%
  merge(.,wGC[,c("qseqid","wX.GC")],by="qseqid",all.x=T)

#order scaffolds
m4<-m4[order(m4$scaffold,m4$start,decreasing = F),]

# tag those windows that are rRNAs
m4$pro.tax[m4$eukRNA.eval>m4$proRNA.eval | (!is.na(m4$proRNA.eval) & is.na(m4$eukRNA.eval))]<-"pro.rRNA"
m4$pro.tax[m4$eukRNA.eval<m4$proRNA.eval | (is.na(m4$proRNA.eval) & !is.na(m4$eukRNA.eval))]<-"euk.rRNA"

# (blastn) Filter all windows with better hit against euk
print("Classifying sliding windows...")

eukWindows<-subset(m4,
                   (euk.bitscore>pro.bitscore  |   # is euk.bitscore higher than pro.bitscore
                    eukRNA.eval<proRNA.eval )  |   # is rRNA evalue higher agains pro than in euk
                   (!is.na(euk.bitscore) & is.na(pro.bitscore)) | # do we have no prok.bitscore but a euk.bitscore?
                   (!is.na(eukRNA.eval) & is.na(proRNA.eval)))    # do we have no prok.rRNA hit but a euk.rRNA hit?

# (blastx) Filter all windows with better hit against euk
eukWindows.x<-subset(m4,
                     euk.x.bitscore>pro.x.bitscore |   # is euk.x.bitscore higher than pro.x.bitscore
                    (!is.na(euk.x.bitscore) & is.na(pro.x.bitscore))) # do we have no prok.x.bitscore but a euk.x.bitscore?

# (blastn) Filter all windows with better hit against pro
proWindows<-subset(m4,   
                   (euk.bitscore<=pro.bitscore  |   # is euk.bitscore lower than pro.bitscore
                    eukRNA.eval>=proRNA.eval )  |   # is rRNA evalue lower agains pro than in euk
                   (is.na(euk.bitscore) & !is.na(pro.bitscore)) | # do we have a prok.bitscore but no euk.bitscore?
                   (is.na(eukRNA.eval) & !is.na(proRNA.eval)))    # do we have a prok.rRNA hit but no euk.rRNA hit?

# (blastx) Filter all windows with better hit against pro
proWindows.x<-subset(m4,
                     euk.x.bitscore<=pro.x.bitscore |   # is euk.x.bitscore lower than pro.x.bitscore
                    (is.na(euk.x.bitscore) & !is.na(pro.x.bitscore))) # do we have a prok.bitscore but no euk.bitscore?

### Identify rate of "prokaryotic windows" across entire chromosom
print("Gathering sliding-window inforation into scafolds...")

# count windows by scaffold
windows<-windowCount %>% 
  dplyr::group_by(V1) %>%
  dplyr::summarise(no_rows = length(V1))
colnames(windows)<-c("scaffold","windowCount")

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

chrSummary<-merge(GC,proWindowCount,by.x="scaffold",by.y="scaffold",all.y=T,all.x=T) %>% 
  merge(.,proWindowCount.x,by.x="scaffold",by.y="scaffold",all.y=T,all.x=T) %>%
  merge(.,eukWindowCount,by.x="scaffold",by.y="scaffold",all.x=T,all.y=T) %>%
  merge(.,eukWindowCount.x,by.x="scaffold",by.y="scaffold",all.x=T,all.y=T)

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

### Retrieve most frequently hit euk_taxon for each scaffold
eukWindowsList<-split(eukWindows,f = eukWindows$scaffold,drop = T)

euk.tmp<-lapply(eukWindowsList, FUN= function(x) {  
  x %>% 
    dplyr::group_by(euk.tax) %>%
    dplyr::summarise(count = length(euk.tax)) %>%
    dplyr::arrange(desc(count)) %>%
    dplyr::slice(1)
})

euk.bestMatch<-plyr::ldply(euk.tmp, rbind)

#calculate coverage across scaffolds
print("Calculating coverage of scaffolds...")

cov$relCov<-log((cov$cov+1)/(median(cov$cov)),2)
covScf<-cov %>% 
  group_by(V1) %>%
  summarise(no_rows = mean(relCov))
colnames(covScf)<-c("scaffold","coverage")

#Combine all scaffold-wide information 
chrSum<-merge(chrSummary,euk.bestMatch[,1:2],by.x="scaffold",by.y=".id",all.x=T) %>%
  merge(.,pro.bestMatch[,1:2],by.x="scaffold",by.y=".id",all.x=T) %>%
  merge(.,covScf,by="scaffold",all.x=T)

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

print("Classifying scaffolds...")

chrSum.euk <- chrSum %>% filter(ratio == 0)
GC_medianCI.euk <- MedianCI(chrSum.euk$GC, method = "boot", conf.level = 0.95)
Coverage_medianCI.euk <- MedianCI(chrSum.euk$coverage, method = "boot", conf.level = 0.95)

chrSum$kingdom<-ifelse(is.na(chrSum$ratio) & is.na(chrSum$bsratio.x), "unknown", 
                ifelse(chrSum$ratio==1, "pro",
                ifelse(chrSum$ratio<=0.5, "euk",
                ifelse(chrSum$bsratio>200, "pro",
                ifelse(chrSum$bsratio>0 & chrSum$bsratio.x>=50 & chrSum$GC<GC_medianCI.euk[2] & chrSum$GC>GC_medianCI.euk[3] &
                       chrSum$coverage<Coverage_medianCI.euk[2] & chrSum$coverage>Coverage_medianCI.euk[3],"pro","uncertain")))))

chrSum$type<-ifelse(chrSum$kingdom=="unknown", "unknown", 
             ifelse(chrSum$kingdom=="uncertain", "uncertain", 
             ifelse(chrSum$kingdom=="pro", chrSum$pro.tax, 
             ifelse(chrSum$kingdom=="euk", chrSum$euk.tax, "ERROR"))))

chrSum$type[is.na(chrSum$type)] <- "unknown"

chrSum$kingdom.bitscore<-ifelse(chrSum$kingdom=="pro", chrSum$pro.bitscore, 
                         ifelse(chrSum$kingdom=="euk", chrSum$euk.bitscore, 0))

# flag scaffolds where the cumulative bitscore is below 200 ("weak evidence")
chrSum$Evidence<-ifelse(chrSum$kingdom.bitscore<200 | chrSum$type=="unknown" | chrSum$type=="uncertain", "weak", "strong")

# select all scaffolds that are identified as prokaryotic
print("Calculating coverage of scaffolds...")

contaminants<-subset(chrSum,kingdom=="pro")
contaminantsTaxaList<-split(contaminants,f=contaminants$pro.tax)
contaminantScaffoldSummary<-contaminants[, c("scaffold", "Length", "GC", "ratio", "bsratio", "ratio.x", "bsratio.x", 
                                                      "pro.tax", "Evidence")]
  
# generate an overview of contaminant scaffolds
contaminantsSummary<-contaminants %>% 
  summarise(Length.mean = round(mean(Length), 3), Length.median = median(Length),
            GC.mean = round(mean(GC),3), GC.median = median(GC))

# generate a simple contaminant table (each raw represents a bac taxon)
contaminantsTable<-contaminants %>% 
                   group_by(pro.tax) %>%
                   summarise(coverage = round(mean(coverage),3),gc=round(mean(GC),3),scaffolds=n())
contaminantsTable$prevalence <- with(contaminantsTable, scaffolds/nrow(chrSum))

# save/print Table for the identified contaminants
print("Generating bacterial scaffolds list/report...")
write.table(contaminantsTable, file = paste(folder,"contaminantsTable.csv",sep=""), quote = F, sep = "\t")
write.table(contaminantScaffoldSummary, file = paste(folder,"contaminantscaffolds.csv",sep=""), quote = F, sep = "\t")

### extract top 10 taxa (to simplify the plot legend)
type.df <- as.data.frame(table(as.factor(chrSum$type)))
type.df.ordered <- type.df[order(-type.df$Freq),]
top10_type.df <- type.df.ordered[1:10,]
top10_type.char <- as.character(top10_type.df$Var1)
top10_type.char <- append(top10_type.char, "other sp.")
chrSum <- chrSum %>% mutate(plot_type = top10_type.char[match(type, top10_type.char, nomatch = 11)])

### Plot pro window percentage vs GC
## kingdom (euk-or-pro) plot
print("Generating scaffold plot...")
pmain<-ggplot(chrSum, aes(x=GC, y=coverage)) + 
  geom_point(aes(size=Length/1e+6,fill=kingdom,color=Evidence),pch=21,alpha=ifelse(chrSum$Evidence=="weak",0.1,.7),stroke=.7)+
  theme_classic()+
  ylab(label = "relative Coverage (log2)")+
  scale_color_manual(name="Evidence",values=c("black","white"))+
  scale_size(name="Size (Mbp)")+
  scale_fill_discrete(name="")+
  theme(legend.spacing = unit(.05,"mm"))
ggsave(paste(folder,"Taxa_screen.kingdom.pdf",sep=""), pmain, device = "pdf")

## top 10 taxa plot
pmain.top10<-ggplot(chrSum, aes(x=GC, y=coverage)) + 
  geom_point(aes(size=Length/1e+6,fill=plot_type,color=Evidence),pch=21,alpha=ifelse(chrSum$Evidence=="weak",0.1,.7),stroke=.7)+
  theme_classic()+
  ylab(label = "relative Coverage (log2)")+
  scale_color_manual(name="Evidence",values=c("black","white"))+
  scale_size(name="Size (Mbp)")+
  scale_fill_discrete(name="")+
  theme(legend.spacing = unit(.05,"mm"))
ggsave(paste(folder,"Taxa_screen.top10taxa.pdf",sep=""), pmain.top10, device = "pdf")

### PCA plot
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
