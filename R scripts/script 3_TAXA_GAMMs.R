rm(list=ls())
'%!in%'<-function(x,y)!('%in%'(x,y))
library(phyloseq)
library(tidyverse)
library(gridExtra)
library(plyr) 
library(reshape2)
library(lubridate)
library(grid)
library(ggplot2)
library(lattice)
library(ggpubr)
library(mgcv)  
library(mgcViz)
library(MASS)
theme_set (theme_classic(base_size=15))


setwd("~/github/")

############################# DATASET #######################################
load("infant and female gelada microbiomes_SILVA_ASV filter 500 reads.RData")
gelada_physeq #3877 ASVs - 1145 samples
otu_table(gelada_physeq)[1:5,1:5]#rows=ASV,columns=samples
head(sample_data(gelada_physeq))
summary(sample_data(gelada_physeq))

#Select only infant samples
table(sample_data(gelada_physeq)$SampleType)
gelada_physeq=prune_samples(sample_data(gelada_physeq)$SampleType %in% "Infant",gelada_physeq)
gelada_physeq=prune_taxa(taxa_sums(gelada_physeq) > 0,gelada_physeq)
gelada_physeq #3784 ASV - 525 samples

meta<-data.frame(sample_data(gelada_physeq))
head(meta)
dim(meta)


### CLEAN SILVA TAXONOMY ###
tax_table(gelada_physeq)<-apply(tax_table(gelada_physeq),2, function(x) (gsub("Unassigned*",NA,x)))
tax_table(gelada_physeq)<-apply(tax_table(gelada_physeq),2, function(x) (gsub("uncultured*",NA, x)))
tax_table(gelada_physeq)<-apply(tax_table(gelada_physeq),2, function(x) (gsub("unidentified*",NA, x)))
tax_table(gelada_physeq)<-apply(tax_table(gelada_physeq),2, function(x) (gsub("gut metagenome",NA, x)))
tax_table(gelada_physeq)<-apply(tax_table(gelada_physeq),2, function(x) (gsub("metagenome",NA, x)))
tax_table(gelada_physeq)<-apply(tax_table(gelada_physeq),2, function(x) (gsub("Mitochondria",NA, x)))
tax_table(gelada_physeq)[,"Family"] [tax_table(gelada_physeq)[,"Order"] %in% "WCHB1-41" ]<- "RFP12" 


### GAMMs per taxa  ###
gamm_per_taxon<- function (dataset,level) {
  
  # Agglomeration at the given taxonomic level
  ps1 = tax_glom(dataset,level,NArm=FALSE)
  ps2<-psmelt(ps1)
  ps2$Taxa<-ps2[, which(colnames(ps2)==level)]
  head(ps2)
  
  # Filter low abundant taxa (<0.01% relative abundance)
  ps2$RelativeAbundance<-(as.numeric(as.character(ps2$Abundance))*100)/as.numeric(as.character(ps2$NumberReads))
  a1<-sapply(unique(ps2$Taxa),function(x) mean(ps2$RelativeAbundance[ps2$Taxa==x],na.rm=TRUE))
  a2<-data.frame(cbind(Taxa=as.character(unique(ps2$Taxa)),RA=a1))
  a2$RA<-as.numeric(as.character(a2$RA))
  ps3<-droplevels(subset(ps2,ps2$Taxa %!in% NA))
  ps3<-droplevels(subset(ps3,ps3$Taxa %in% unique(a2$Taxa[a2$RA>=0.01])))

  # Loop
  df2<-NULL
  for (i in unique(ps3$Taxa)){
    print(i)
    data<-droplevels(subset(ps3,ps3$Taxa %in% i))
    m1<-gam(log10(RelativeAbundance+0.001) ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30 , data=data, gamma=1.4,method="REML") 
    coef<-summary(m1)$p.table[,1]
    coef1<-summary(m1)$s.table[,1]
    pval<-summary(m1)$p.table[,4]
    pval1<-summary(m1)$s.table[,4]
    dev<-summary(m1)$dev.expl*100
    df<-c(i,dev,coef,coef1,pval,pval1)
    df2<-rbind(df2,df)
    colnames(df2)<-c("Taxa","Deviance","coefIntercept","coefSexM","coefParityPrimi","coefMomRank","coefRain30","coefMinT30","s.Age","s.Individual","s.Unit","pIntercept","pSex","pParity","pMomRank","pRain30","pMinT30","pAge","pIndividual","pUnit")
  }
  df2<-data.frame(df2)
  
  #Adjust p-values for multiple testing
  df2$pIntercept<-as.numeric(as.character(df2$pIntercept))
  df2$pSex<-as.numeric(as.character(df2$pSex))
  df2$pParity<-as.numeric(as.character(df2$pParity))
  df2$pMomRank<-as.numeric(as.character(df2$pMomRank))
  df2$pRain30<-as.numeric(as.character(df2$pRain30))
  df2$pMinT30<-as.numeric(as.character(df2$pMinT30))
  df2$pAge<-as.numeric(as.character(df2$pAge))
  df2$pIndividual<-as.numeric(as.character(df2$pIndividual))
  df2$pUnit<-as.numeric(as.character(df2$pUnit))
  
  df2$padjIntercept<-p.adjust(df2$pIntercept,method="fdr")
  df2$padjSex<-p.adjust(df2$pSex,method="fdr")
  df2$padjParity<-p.adjust(df2$pParity,method="fdr")
  df2$padjMomRank<-p.adjust(df2$pMomRank,method="fdr")
  df2$padjRain30<-p.adjust(df2$pRain30,method="fdr")
  df2$padjMinT30<-p.adjust(df2$pMinT30,method="fdr")
  df2$padjAge<-p.adjust(df2$pAge,method="fdr")
  df2$padjIndividual<-p.adjust(df2$pIndividual,method="fdr")
  df2$padjUnit<-p.adjust(df2$pUnit,method="fdr")
  
  return(df2)
}


### ALL SAMPLES ###
res<-gamm_per_taxon(gelada_physeq,"Family")
#res<-gamm_per_taxon(gelada_physeq,"Genus")

### ONLY ON YOUNG IMMATURES ###
young=prune_samples(sample_data(gelada_physeq)$AgeMonth<=12, gelada_physeq)
young=prune_taxa(taxa_sums(young) > 0,young)
young #184 samples
#res<-gamm_per_taxon(young,"Family")
#res<-gamm_per_taxon(young,"Genus")

### ONLY ON OLD IMMATURES ###
old=prune_samples(sample_data(gelada_physeq)$AgeMonth>18, gelada_physeq)
old=prune_taxa(taxa_sums(old) > 0,old)
old #259 samples
#res<-gamm_per_taxon(old,"Family")
#res<-gamm_per_taxon(old,"Genus")



### CHECK RES ##
head(res)
dim(res[res$padjAge<=0.05,])
dim(res[res$padjParity<=0.05,])
dim(res[res$pParity<=0.05,])#non-adjusted p-values
res[res$pParity<=0.05,]
dim(res[res$padjMomRank<=0.05,])
dim(res[res$padjGrSize<=0.05,])
dim(res[res$padjSex<=0.05,])
dim(res[res$padjRain30<=0.05,])
dim(res[res$padjMinT30<=0.05,])

#Plot specific taxa
ps1 = tax_glom(gelada_physeq,"Family",NArm=FALSE) #CHANGE LEVEL HERE
ps2<-psmelt(ps1)
ps2$RelativeAbundance<-(as.numeric(as.character(ps2$Abundance))*100)/as.numeric(as.character(ps2$NumberReads))
ps2$Taxa<-ps2$Family #CHANGE LEVEL HERE
head(ps2)
data<-droplevels(subset(ps2,ps2$Taxa %in% "Desulfovibrionaceae"))
ggplot(data, aes(x=AgeMonth, y=RelativeAbundance+0.001,colour=Parity)) + geom_point()+scale_y_log10() + ylab("Log relative abundance (%)")+ xlab("Age (months)")+geom_smooth(method="gam",se=F)




