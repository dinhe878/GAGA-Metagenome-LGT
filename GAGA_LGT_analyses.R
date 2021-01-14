#!/usr/bin/env Rscript

## Load packages
library(cowplot)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(GenomicRanges)
library(ggfortify)
library(ggpubr)
library(DescTools)

## Prepare files

# GAGA ID
args <- commandArgs(trailingOnly = T)
if (length(args)==0) { args[1] <- "GAGA-0024" }
id<-args[1]

# local result folder
folder<-paste("/Users/dinghe/Dropbox/Work/Projects/GAGA/Metagenome_analysis/LGT_pipeline/Results/Latest/14012021/",id,"/",sep="")

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

## LGT Analysis
# blastn vs prok database (NCBI: 1908 complete prokaryotic genomes)
# blastn vs insect database (NCBI: 43 high quality complete genomes)
# blastx vs prok database (Uniprot90 database - Bacteria part)
# blastx vs insect database (Uniprot90 database - Insecta part)
# rRNA annotations
# genomic coverage from long-read data
# GC content 


# name columns
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

#euk$euk.tax[euk$euk.tax=="unclassified"]<-"Insecta"

euk.x$euk.x.tax <- with(euk.x, ifelse(grepl(".*g_(.*?);.*", euk.x.sgi), gsub(".*g_(.*?);.*","\\1",euk.x.sgi), 
                             ifelse(grepl(".*f_(.*?);.*", euk.x.sgi), gsub(".*f_(.*?);.*","\\1",euk.x.sgi),
                             ifelse(grepl(".*o_(.*?);.*", euk.x.sgi), gsub(".*o_(.*?);.*","\\1",euk.x.sgi),
                             ifelse(grepl(".*c_(.*?);.*", euk.x.sgi), gsub(".*c_(.*?);.*","\\1",euk.x.sgi),
                             ifelse(grepl(".*p_(.*?);.*", euk.x.sgi), gsub(".*p_(.*?);.*","\\1",euk.x.sgi),
                                                                    gsub(".*d_(.*?);.*","\\1",euk.x.sgi)))))))
#euk.x$euk.x.tax[euk.x$euk.x.tax=="unclassified"]<-"Insecta"

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

# Compute alignment lenghts
# Currently not in use
#m4$eukAln<-abs(m4$euk.send-m4$euk.sstart)
#m4$proAln<-abs(m4$pro.send-m4$pro.sstart)

# (blastn) Filter all windows with better hit against euk
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

### Identify rate of "prokaryotic windows" across entire chromosome
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
# we likely do not need to merge "windows" df?
#%>% merge(.,windows,by.x="scaffold",by.y="scaffold",all.x=T,all.y=T) 

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
cov$relCov<-log((cov$cov+1)/(median(cov$cov)),2)
covScf<-cov %>% 
  group_by(V1) %>%
  summarise(no_rows = mean(relCov))
colnames(covScf)<-c("scaffold","coverage")

#Combine all scaffold-wide information 
chrSum<-merge(chrSummary,euk.bestMatch[,1:2],by.x="scaffold",by.y=".id",all.x=T) %>%
  merge(.,pro.bestMatch[,1:2],by.x="scaffold",by.y=".id",all.x=T) %>%
  merge(.,covScf,by="scaffold",all.x=T)

# flag a scaffold as prokaryotic if 
#   1. 100 % of hitted windows are identified as prokaryotic OR
#   2. more than 50% of hitted windows are identified as prokaryotic AND bsratio > 200
#   3. more than 50% of hitted windows are identified as prokaryotic AND 0 < bsratio < 200 AND thier
#     3.1 bsratio.x > 50 AND
#     3.2 GC content does not lie within 95% confidence interval of euk scaffolds GC content distribution AND
#     3.3 coverage does not lie within 95% confidence interval of euk scaffolds coverage distribution 
# note that it is possible to have ratio 1 for only small amount of proWindow but the rest are uncertain (no hit)
# for downstream genome annotation, only the scaffolds tagged "pro" should be excluded
# for downstram microbiome analyses, include all "pro" and "uncertain" scaffolds

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
contaminants<-subset(chrSum,kingdom=="pro")
contaminantsList<-split(contaminants,f=contaminants$pro.tax)

# generate an overview of contaminant scaffolds
contaminantsSummary<-contaminants %>% 
  summarise(Length.mean = round(mean(Length), 3), Length.median = median(Length),
            GC.mean = round(mean(GC),3), GC.median = median(GC), 
            )

# generate a simple contaminant table (each raw represents a bac taxon)
contaminantsTable<-contaminants %>% 
  group_by(pro.tax) %>%
  summarise(coverage = round(mean(coverage),3),gc=round(mean(GC),3),scaffolds=n())

# save/print Table for the identified contaminants
write.table(contaminantsTable, file = paste(folder,"contaminantsTable.csv",sep=""), quote = F, sep = "\t")
#contaminantsTable

### extract top 10 taxa (to simplify the plot legend)
type.df <- as.data.frame(table(as.factor(chrSum$type)))
type.df.ordered <- type.df[order(-type.df$Freq),]
top10_type.df <- type.df.ordered[1:10,]
top10_type.char <- as.character(top10_type.df$Var1)
top10_type.char <- append(top10_type.char, "other sp.")
chrSum <- chrSum %>% mutate(plot_type = top10_type.char[match(type, top10_type.char, nomatch = 11)])

### Plot pro window percentage vs GC
## kingdom (euk-or-pro) plot
pmain<-ggplot(chrSum, aes(x=GC, y=coverage)) + 
  geom_point(aes(size=Length/1e+6,fill=kingdom,color=Evidence),pch=21,alpha=ifelse(chrSum$Evidence=="weak",0.1,.7),stroke=.7)+
  theme_classic()+
  ylab(label = "relative Coverage (log2)")+
  scale_color_manual(name="Evidence",values=c("black","white"))+
  scale_size(name="Size (Mbp)")+
  scale_fill_discrete(name="")+
  theme(legend.spacing = unit(.05,"mm"))
legend<-get_legend(pmain)
pmain<-pmain+theme(legend.position = "none")

# Main plot
# Marginal densities along x axis
xstats <- chrSum %>% group_by(type) %>% summarise(median = median(GC),n = n())
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = chrSum, aes(x = GC, fill = type),
               alpha = 0.7, size = 0.2)
#+geom_vline(data=xstats,aes(xintercept=median,color=type),size=.2,
#              linetype="dashed")+scale_size_continuous(range=c(0,1))

# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ystats <- chrSum %>% group_by(type) %>% summarise(median = median(coverage),n = n())
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = chrSum, aes(x = coverage, fill = type),
               alpha = 0.7, size = 0.2) + coord_flip()
#+geom_vline(data=ystats,aes(xintercept=median,color=type),size=.2,
#              linetype="dashed")+scale_size_continuous(range=c(0,1))

p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3<-ggdraw(p2)

#tab<-ggtexttable(contaminantsTable, rows = NULL, 
#                        theme = ttheme("classic",base_size = 7,padding = unit(c(1, 1), "mm")))

p4 <- plot_grid(p3,legend,nrow=2,ncol=2,rel_widths = c(.75,.25),rel_heights = c(0.8,.2))
ggsave(paste(folder,"LGT_screen.kingdom.pdf",sep=""), p4, device = "pdf")
#dev.print(pdf,paste(folder,"proScreen.pdf",sep=""),width=10,height=7)


## prettyer top 10 taxa plot
pmain.top10<-ggplot(chrSum, aes(x=GC, y=coverage)) + 
  geom_point(aes(size=Length/1e+6,fill=plot_type,color=Evidence),pch=21,alpha=ifelse(chrSum$Evidence=="weak",0.1,.7),stroke=.7)+
  theme_classic()+
  ylab(label = "relative Coverage (log2)")+
  scale_color_manual(name="Evidence",values=c("black","white"))+
  scale_size(name="Size (Mbp)")+
  scale_fill_discrete(name="")+
  theme(legend.spacing = unit(.05,"mm"))
legend.top10<-get_legend(pmain.top10)
pmain.top10<-pmain.top10+theme(legend.position = "none")
p1 <- insert_xaxis_grob(pmain.top10, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3<-ggdraw(p2)
p4 <- plot_grid(p3,legend.top10,nrow=2,ncol=2,rel_widths = c(.75,.25),rel_heights = c(0.8,.2))
ggsave(paste(folder,"LGT_screen.top10taxa.pdf",sep=""), p4, device = "pdf")






### To-do $ experiments
# tag NA to pro-like or euk-like
# for scaffolds not hitting anything -> tag as ant
# final strategy for taging bac scaffolds
#   Criteria to assign a bac scaffold
#   1. ratio == 1
#   2. (0.5 < ratio < 1) & GC clustering does not belongs to the bulk of euk scaffolds (ratio == 0)
#      & Coverage clustering belongs to the bulk of euk scaffolds (ratio == 0) 
#        

### Distribution analyses
library(fitdistrplus)
library(gamlss)
library(gamlss.dist)
library(gamlss.add)
library(DescTools)

ggdensity(chrSum.euk$GC, main = "Euk GC content")
ggdensity(chrSum.pro$GC, main = "Pro GC content")
ggqqplot(chrSum.euk$GC)
ggqqplot(chrSum.pro$GC)
shapiro.test(chrSum.euk$GC)
descdist(chrSum.euk$GC)


fit <- fitDist(chrSum.pro$GC, k = 2, type = "realline", trace = F, try.gamlss = T)
summary(fit)
plot(fit)
histDist(chrSum.pro$GC, "SEP2")
wp(fit)
prof.dev(fit, "mu", min = 1, max = 99)

### PCA analyses
library(factoextra)
chrSum.pca_df <- chrSum[c("scaffold", "GC", "ratio", "coverage", "kingdom", "Evidence")] %>% 
  tibble::column_to_rownames(var = "scaffold") %>%
  na.omit()

# mutate(Evidence = ifelse(Evidence == "strong", 1, 0.1)) %>%
# chrSum.pca_df.ambiguous_cases <- chrSum.pca_df %>% filter(ratio > 0.5 & ratio < 1)

chrSum.pca_df.pca <- dudi.pca(chrSum.pca_df[1:3], nf=3, scannf = FALSE)
p.pca <- fviz_pca_biplot(chrSum.pca_df.pca, label = c("var", "quali"), 
                fill.ind = chrSum.pca_df$kingdom, col.ind = chrSum.pca_df$Evidence, 
                pointshape = 21, palette = c("red","blue","white"),
                title = id, legend.title = list(color = "Evidence", fill = "Taxa"), mean.point = F 
                )
p.pca


# pointsize = chrSum.pca_df$Evidence, 
