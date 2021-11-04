rm(list=ls())
'%!in%'<-function(x,y)!('%in%'(x,y))
library(phyloseq)
library(tidyverse)
library(vegan)
library(gridExtra)
library(reshape2)
library(lubridate)
library(ggplot2)
library(compositions)
theme_set (theme_classic(base_size=15))


setwd("~/github/")

############################# DATASET #######################################
load("datasets/infant and female gelada microbiomes_SILVA_ASV filter.RData")
gelada_physeq #3877 ASV - 1145 samples
otu_table(gelada_physeq)[1:5,1:5]#rows=ASV,columns=samples
head(sample_data(gelada_physeq))
summary(sample_data(gelada_physeq))

#Infant samples with a maternal match
INF = prune_samples(sample_data(gelada_physeq)$SampleType %in% "Infant" & sample_data(gelada_physeq)$MatchedMomSample %!in% NA, gelada_physeq)
INF #398 

#Female samples
AF = prune_samples(sample_data(gelada_physeq)$SampleType %in% "Adult female",gelada_physeq) 
AF #620

dataset<-merge_phyloseq(INF,AF)

#Rarefy datatset
set.seed(123)
raref <- rarefy_even_depth(dataset,sample.size=20000,replace=FALSE,trimOTUs=TRUE)
otu_table(raref)[1:5,1:5]
colSums(otu_table(raref))
meta<-data.frame(sample_data(raref))
head(meta)


################### Calculate beta dissimilarity between all pairs of samples ###################
#Bray-Curtis distance
bcdist = phyloseq::distance(raref, method="bray",normalized=TRUE, parallel=FALSE, fast=TRUE)  
bcdist1=as.matrix(bcdist)
bcdist2<-melt(bcdist1)
colnames(bcdist2)[3]<-"bcval"
head(bcdist2)

#Unweighted UniFrac distance
udist=UniFrac(raref,weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)
udist1=as.matrix(udist)
udist2<-melt(udist1)
colnames(udist2)[3]<-"uval"

#Weighted UniFrac distance
wdist=UniFrac(raref , weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
wdist1=as.matrix(wdist)
wdist2<-melt(wdist1)
colnames(wdist2)[3]<-"wval"

### Add distance to dataset1 ###
bcdist2$Match<-paste(bcdist2$Var1,bcdist2$Var2,sep="")
udist2$Match<-paste(udist2$Var1,udist2$Var2,sep="")
wdist2$Match<-paste(wdist2$Var1,wdist2$Var2,sep="")
beta1<-merge(bcdist2,udist2[c("Match","uval")],by="Match",all.x=TRUE,all.y=FALSE)
beta2<-merge(beta1,wdist2[c("Match","wval")],by="Match",all.x=TRUE,all.y=FALSE)
head(beta2)
dim(beta2)
colnames(beta2)[2:3]<-c("Individual1","Individual2")

#Add metadata of sample 1
beta3<-merge(beta2,meta[c("SampleID","Individual","SampleType","Date","Sex","Unit","AgeMonth","Season")],by.x="Individual1",by.y="SampleID",all.x=TRUE,all.y=FALSE)
head(beta3)

#Add metadata of sample 2
beta4<-merge(beta3,meta[,c("SampleID","Individual","SampleType","Date","Sex","Unit","AgeMonth","Season")],by.x="Individual2",by.y="SampleID",all.x=TRUE,all.y=FALSE)
head(beta4)

beta4$Same<-ifelse(as.character(beta4$Individual1)==as.character(beta4$Individual2),1,0)
beta5<-droplevels(subset(beta4,beta4$Same %!in% 1))
beta6<-droplevels(subset(beta5,beta5$SampleType.x %in% "Infant" & beta5$SampleType.y %in% "Adult female"))
head(beta6)
dim(beta6)

# Number of shared ASVs
nb1<-NULL
Shared<-list()
Name<-list()
for (j in 1:length(beta6$Match)){
  print(beta6$Match[j])
  #List ASVs from the first sample
  col1<-which(colnames(otu_table(raref)) %in% beta6$Individual1[j])
  s1<-data.frame(otu_table(raref)[,col1])
  colnames(s1)<-"Sample1"
  s1<-droplevels(subset(s1,s1$Sample1 %!in% 0))
  s1.NumberASV<-nrow(s1)
  s1.list<-rownames(s1)
  
  #List ASVs from the second sample
  col2<-which(colnames(otu_table(raref)) %in% beta6$Individual2[j])
  s2<-data.frame(otu_table(raref)[,col2])
  colnames(s2)<-"Sample2"
  s2<-droplevels(subset(s2,s2$Sample2 %!in% 0))
  s2.NumberASV<-nrow(s2)
  s2.list<-rownames(s2)
  NumberSharedASV<-length(intersect(s1.list,s2.list))
  nb<-cbind(as.character(beta6$Match[j]),s1.NumberASV,s2.NumberASV,NumberSharedASV)
  nb1<-rbind(nb1,nb)
  Shared[[j]]<-intersect(s1.list,s2.list)
  Name[[j]]<-rep(as.character(beta6$Match[j]),length(Shared[[j]]))
}
nb1<-data.frame(nb1)
head(nb1)
dim(nb1)

beta7<-merge(beta6,nb1,by.y="V1",by.x="Match",all.x=TRUE,all.y=TRUE)
head(beta7)

beta7$s1.NumberASV <-as.numeric(as.character(beta7$s1.NumberASV))
beta7$s2.NumberASV <-as.numeric(as.character(beta7$s2.NumberASV))
beta7$NumberSharedASV <-as.numeric(as.character(beta7$NumberSharedASV))
beta7$PropSharedASV<-(beta7$NumberSharedASV*100)/beta7$s1.NumberASV
summary(beta7)
names(beta7)
colnames(beta7)[c(2,3,7:20,22:23)]<-c("SampleID2","SampleID1",
                                      "Ind1","SampleType1","Date1","Sex1","Unit1","AgeMonth1","Season1",
                                      "Ind2","SampleType2","Date2","Sex2","Unit2","AgeMonth2","Season2",
                                       "NumberASV1","NumberASV2")

beta7$SameSeason<-ifelse(as.character(beta7$Season1)==as.character(beta7$Season2),1,0)   
beta7$SameUnit<-ifelse(as.character(beta7$Unit1)==as.character(beta7$Unit2),1,0)   
beta7$Match<-as.factor(beta7$Match)
beta7$Date1<-as.Date(beta7$Date1,"%Y-%m-%d")
beta7$Date2<-as.Date(beta7$Date2,"%Y-%m-%d")
beta7$SameSeason<-as.factor(beta7$SameSeason)
beta7$SameUnit<-as.factor(beta7$SameUnit)
summary(beta7)
names(beta7)
head(beta7)
head(beta7[c(1,3,2,7:20,26:27,4:6,22:25)])
#saveRDS(beta7[c(1,3,2,7:20,26:27,4:6,22:25)],"datasets/mother-offspring pairs/CompositionalSimilarity_ALL PAIRS.txt")

#Taxa shared between samples
TaxaShared<-data.frame(cbind(unlist(Name),unlist(Shared)))
head(TaxaShared)
dim(TaxaShared)
colnames(TaxaShared)<-c("Match","ASVShared")

#Select only actual mother-offspring pairs
pair<-readRDS("datasets/metadata.txt")
pair<-droplevels(subset(pair,pair$MatchedMomSample %!in% NA))
pair$Match<-paste(pair$SampleID,pair$MatchedMomSample,sep="")
head(pair)
dim(pair) #398 mother-infant pairs (collected0 or 1 day apart) 
TaxaShared1<-droplevels(subset(TaxaShared,TaxaShared$Match %in% pair$Match))
dim(TaxaShared1)

#Add metadata
TaxaShared2<-merge(TaxaShared1,beta7,by="Match",all.x=TRUE,all.y=FALSE)
head(TaxaShared2)
dim(TaxaShared2)
summary(TaxaShared2)
length(unique(TaxaShared2$Match))#398 pairs
summary(as.numeric(table(TaxaShared2$Match)))
#saveRDS(TaxaShared2[,-22],"datasets/mother-offspring pairs/TaxaShared_ALL PAIRS.txt")


