rm(list=ls())
'%!in%'<-function(x,y)!('%in%'(x,y))
library(reshape)
library(vegan)
library(ape)
library(bcRep)
library(mgcv)  
library(rgl)
library(mgcViz)
library(MASS)
library(ggplot2)
library(ggpubr)
library(gridExtra)
theme_set (theme_classic(base_size=15))

setwd("~/github/")

############################### PICRUSt2 DATASETS ######################################################
# KEGG LEVEL 2 #
ko2<-read.csv("datasets/picrust2/ko_level2_relative abundance.txt",header=TRUE, row.names=1, sep=",", na.strings="NA") 
dim(ko2)#39
ko2[1:5,1:5]

# KEGG LEVEL 3 #
ko3<-read.csv("datasets/picrust2/ko_level3_relative abundance.txt",header=TRUE, row.names=1,sep=",", na.strings="NA") 
dim(ko3)#271
ko3[1:5,1:5]

# ENZYME NUMBER #
ec<-read.csv("datasets/picrust2/EC_relative abundance.txt",header=TRUE, row.names=1, sep=",", na.strings="NA") 
ec[1:5,1:5]
dim(ec)#1931


### GAMMs per pathway  ###
metadata<-readRDS("datasets/metadata.txt")
metadata<-droplevels(subset(metadata,metadata$SampleType %in% "Infant"))
head(metadata)

gamm_per_pathway<- function (dataset,min.age,max.age) {
  meta<-droplevels(subset(metadata,metadata$AgeMonth>=min.age & metadata$AgeMonth<=max.age))
  dataset<-droplevels(dataset[,which(colnames(dataset) %in%  meta$SampleID)])
  dataset<-droplevels(dataset[which(rowMeans(dataset)>=0.1),])
  dataset$Pathway<-rownames(dataset)
  dataset1<-melt(dataset)
  colnames(dataset1)<-c("Pathway","SampleID","RA")
  dataset2<-merge(dataset1,meta,by="SampleID",all.x=TRUE,all.y=FALSE)

df2<-NULL
for (i in unique(dataset2$Pathway)){
  print(i)
  data<-droplevels(subset(dataset2,dataset2$Pathway %in% i))
  m1<-gam(RA ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30, data=data,gamma=1.4,method="REML") 
  coef<-summary(m1)$p.table[,1]
  coef1<-summary(m1)$s.table[,1]
  pval<-summary(m1)$p.table[,4]
  pval1<-summary(m1)$s.table[,4]
  dev<-summary(m1)$dev.expl*100
  df1<-c(i,dev,coef,coef1,pval,pval1)
  df2<-rbind(df2,df1)
  names(df2)
  colnames(df2)[c(1:2,12:20)]<-c("Pathway","Deviance","pIntercept","pSex","pParity","pMomRank","pRain","pMinT","pAge","pIndividual","pUnit")
}
df2<-data.frame(df2)

#Adjust pvalues for multiple testing
df2$pIntercept<-as.numeric(as.character(df2$pIntercept))
df2$pSex<-as.numeric(as.character(df2$pSex))
df2$pParity<-as.numeric(as.character(df2$pParity))
df2$pMomRank<-as.numeric(as.character(df2$pMomRank))
df2$pRain<-as.numeric(as.character(df2$pRain))
df2$pMinT<-as.numeric(as.character(df2$pMinT))
df2$pAge<-as.numeric(as.character(df2$pAge))
df2$pIndividual<-as.numeric(as.character(df2$pIndividual))
df2$pUnit<-as.numeric(as.character(df2$pUnit))

df2$padjIntercept<-p.adjust(df2$pIntercept,method="fdr")
df2$padjSex<-p.adjust(df2$pSex,method="fdr")
df2$padjParity<-p.adjust(df2$pParity,method="fdr")
df2$padjMomRank<-p.adjust(df2$pMomRank,method="fdr")
df2$padjRain<-p.adjust(df2$pRain,method="fdr")
df2$padjMinT<-p.adjust(df2$pMinT,method="fdr")
df2$padjAge<-p.adjust(df2$pAge,method="fdr")
df2$padjIndividual<-p.adjust(df2$pIndividual,method="fdr")
df2$padjUnit<-p.adjust(df2$pUnit,method="fdr")

return(df2)

} 


### ALL IMMATURES ###
res<-gamm_per_pathway(ko2,0,36)
#res<-gamm_per_pathway(ko3,0,36)
#res<-gamm_per_pathway(ec,0,36)
head(res)
dim(res[res$padjParity<=0.05,])

### ONLY YOUNG IMMATURES ###
#res<-gamm_per_pathway(ko2,0,12)
#res<-gamm_per_pathway(ko3,0,12) 
#res<-gamm_per_pathway(ec,0,12)

### ONLY OLD IMMATURES ###
#res<-gamm_per_pathway(ko2,18,36)
#res<-gamm_per_pathway(ko3,18,36) 
#res<-gamm_per_pathway(ec,18,36)

