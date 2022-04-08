rm(list=ls())
'%!in%'<-function(x,y)!('%in%'(x,y))
library(ggplot2)
library(mgcv)  
library(mgcViz)
library(lubridate)
theme_set (theme_classic(base_size=15))

setwd("~/github/")

###### Dataset of compositional similarity (or dissimilarity) between all possible pairs of infant-adult female samples ######
diss<-readRDS("datasets/mother-offspring pairs/CompositionalSimilarity_ALL PAIRS_css_ASV 500 reads.txt")
head(diss) 
dim(diss)#325.488
summary(diss)
diss$Ind1<-as.factor(diss$Ind1)
diss$Ind2<-as.factor(diss$Ind2)
table(diss$SampleType1) #SampleID1 = always an infant sample
table(diss$SampleType2) #SampleID2 = always an adult female sample
diss$TotReads<-diss$NumberReads1+diss$NumberReads2


######################### OBSERVED VALUE (actual mother-offspring matched pairs) ######################### 
pair<-readRDS("datasets/metadata.txt")
pair<-droplevels(subset(pair,pair$MatchedMomSample %!in% NA))
pair$Match<-paste(pair$SampleID,pair$MatchedMomSample,sep="")
head(pair)
dim(pair)#398 mother-infant pairs (collected 0 or 1 day apart) 

#Add dissimilarity
head(diss)
observed<-merge(pair,diss[c("Match","SameUnit","ReproState2","NumberReads2","DiffRead","DiffDay","bcval","uval","wval","NumberASV1","NumberASV2","NumberSharedASV","PropSharedASV.inf","PropSharedASV.all")],by="Match",all.x=TRUE,all.y=FALSE)
head(observed)
colnames(observed)[27:28]<-c("ReproStateMom","NumberReadsMom")
table(observed$SampleType)


#Only take infants for which we have at least one matched mother-offspring sample
diss<-droplevels(subset(diss,diss$SampleID1 %in% observed$SampleID))
dim(diss)#246.760        



######################### RANDOM POOL ################################# 
#Add info of the actual mother matched sample
diss1<-merge(diss,observed[c("SampleID","Mother","MatchedMomSample","ReproStateMom")],by.x="SampleID1",by.y="SampleID",all.x=TRUE,all.y=FALSE)
head(diss1)
dim(diss1)
diss1$SameRepro<-as.factor(ifelse(as.character(diss1$ReproState2)==as.character(diss1$ReproStateMom),1,0))
table(diss1$SameRepro)

#Add info of observed pairs and maternal samples
diss1$ObservedPair<-ifelse(as.character(diss1$SampleID2)==as.character(diss1$MatchedMomSample),"1","0")
table(diss1$ObservedPair)
diss1$ObservedPair<-as.factor(diss1$ObservedPair)
diss1$MotherSample<-ifelse(as.character(diss1$Ind2)==as.character(diss1$Mother),"1","0")
table(diss1$MotherSample)

#Remove any maternal samples from random pool 
pool0<-droplevels(subset(diss1,diss1$MotherSample %!in% "1")) #remove all matched mother samples from the pool of random females
dim(pool0)#240.429

#Only samples within 60 days of each other
pool1<-droplevels(subset(pool0,pool0$DiffDay<=60))
dim(pool1)

#Only from same unit 
pool2<-droplevels(subset(pool1,pool1$SameUnit %in% 1))
dim(pool2)

#In same repro state than mother
pool3<-droplevels(subset(pool2,pool2$SameRepro %in% 1))
dim(pool3)#3112

length(unique(pool3$SampleID1))#345 pairs left
summary(as.numeric(table(pool3$SampleID1)))#7-9 female samples per infant
length(table(pool3$SampleID1)[table(pool3$SampleID1)<=1])#22
hist(table(pool3$SampleID1),nclass=50,xlab="# samples per infant")


#Set initial observed datasets 
observed<-droplevels(subset(diss1,diss1$ObservedPair %in% "1")) 
dim(observed)#398

observed_young<-droplevels(subset(observed,observed$AgeMonth1<=12))
dim(observed_young)#136

observed_old<-droplevels(subset(observed,observed$AgeMonth1>=18))
dim(observed_old)#201


#### RESAMPLING ####
loop<-function(pool_rand,n_simu){
  
  #Restrict to infant samples with enough resolution 
  observed1<-droplevels(subset(observed,observed$SampleID1 %in% unique(pool_rand$SampleID1))) 
  observed_young1<-droplevels(subset(observed1,observed1$AgeMonth1<=12))
  observed_old1<-droplevels(subset(observed1,observed1$AgeMonth1>=18))
  
  all<-NULL
  young<-NULL
  old<-NULL
  old_sub1<-NULL
  
  for (j in 1:n_simu) {
    print(j)
    
    #Subsample in the old samples to have only same than in young sample
    sub<-sample(observed_old1$SampleID1,size=nrow(observed_young1),replace=FALSE)
    sub_old<-droplevels(subset(observed_old1,observed_old1$SampleID1 %in% sub))
    
    ### Generate random pool of females ###
    rand<-NULL
    rand<-sapply(unique(pool_rand$SampleID1),function(x) sample(pool_rand$SampleID2[pool_rand$SampleID1==x],1,replace=FALSE))
    rand1<-data.frame(cbind(SampleID1=as.character(unique(pool_rand$SampleID1)),SampleID2=as.character(rand)))
    rand1$Match<-paste(rand1$SampleID1,rand1$SampleID2,sep="")
    rand2<-merge(rand1,diss1[,-c(1,3)],by="Match",all.x=TRUE,all.y=FALSE)# SampleID1 + SampleID2 
    dataset<-rbind(observed1,rand2)
    
    all_pair<-dataset
    young_pair<-droplevels(subset(dataset,dataset$AgeMonth1<=12 & dataset$SampleID1 %in% observed_young1$SampleID1))
    old_pair<-droplevels(subset(dataset,dataset$AgeMonth1>=18 & dataset$SampleID1 %in% observed_old1$SampleID1))
    old_pair_sub<-droplevels(subset(dataset,dataset$AgeMonth1>=18 & dataset$SampleID1 %in% sub_old$SampleID1))
    
    mod<-function(dataset1){
      res<-NULL
      obs<-droplevels(subset(dataset1,dataset1$ObservedPair %in% 1))
      oMeanNumberSharedASV<-mean(obs$NumberSharedASV)
      oMeanPropSharedASV.inf<-mean(obs$PropSharedASV.inf)
      oMeanPropSharedASV.all<-mean(obs$PropSharedASV.all)
      oMeanBC<-mean(obs$bcval)
      oMeanUU<-mean(obs$uval)
      oMeanWU<-mean(obs$wval)
      
      random<-droplevels(subset(dataset1,dataset1$ObservedPair %in% 0))
      rMeanNumberSharedASV <-mean(random$NumberSharedASV)
      rMeanPropSharedASV.inf<-mean(random$PropSharedASV.inf)
      rMeanPropSharedASV.all<-mean(random$PropSharedASV.all)
      rMeanBC<-mean(random$bcval)
      rMeanUU<-mean(random$uval)
      rMeanWU<-mean(random$wval)
      
      dataset1$ObservedPair<-factor(dataset1$ObservedPair,levels = c("0","1"))
      m1<-gam(NumberSharedASV ~ ObservedPair + s(AgeMonth1) + Sex1 + TotReads + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset1,gamma=1.4,method="REML") 
      m2<-gam(PropSharedASV.inf ~ ObservedPair + s(AgeMonth1) + Sex1 + TotReads + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset1,gamma=1.4,method="REML") 
      m3<-gam(PropSharedASV.all ~ ObservedPair + s(AgeMonth1) + Sex1 + TotReads + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset1,gamma=1.4,method="REML") 
      m4<-gam(bcval ~ ObservedPair + s(AgeMonth1) + Sex1  + TotReads + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset1,gamma=1.4,method="REML") 
      m5<-gam(uval ~ ObservedPair + s(AgeMonth1) + Sex1 + TotReads + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset1,gamma=1.4,method="REML") 
      m6<-gam(wval ~ ObservedPair+ s(AgeMonth1) + Sex1 + TotReads + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset1,gamma=1.4,method="REML") 
      
      jd1<-c(MeanObs=oMeanNumberSharedASV,MeanRandom=rMeanNumberSharedASV,summary(m1)$p.table[,1],summary(m1)$s.table[,1],summary(m1)$p.table[,4],summary(m1)$s.table[,4],summary(m1)$dev.expl*100,NbDyads=nrow(dataset1)/2,Type="NumberSharedASV")
      jd2<-c(MeanObs=oMeanPropSharedASV.inf,MeanRandom=rMeanPropSharedASV.inf,summary(m2)$p.table[,1],summary(m2)$s.table[,1],summary(m2)$p.table[,4],summary(m2)$s.table[,4],summary(m2)$dev.expl*100,NbDyads=nrow(dataset1)/2,Type="PropSharedASV.inf")
      jd3<-c(MeanObs=oMeanPropSharedASV.all,MeanRandom=rMeanPropSharedASV.all,summary(m3)$p.table[,1],summary(m3)$s.table[,1],summary(m3)$p.table[,4],summary(m3)$s.table[,4],summary(m3)$dev.expl*100,NbDyads=nrow(dataset1)/2,Type="PropSharedASV.all")
      jd4<-c(MeanObs=oMeanBC,MeanRandom=rMeanBC,summary(m4)$p.table[,1],summary(m4)$s.table[,1],summary(m4)$p.table[,4],summary(m4)$s.table[,4],summary(m4)$dev.expl*100,NbDyads=nrow(dataset1)/2,Type="BC")
      jd5<-c(MeanObs=oMeanUU,MeanRandom=rMeanUU,summary(m5)$p.table[,1],summary(m5)$s.table[,1],summary(m5)$p.table[,4],summary(m5)$s.table[,4],summary(m5)$dev.expl*100,NbDyads=nrow(dataset1)/2,Type="UU")
      jd6<-c(MeanObs=oMeanWU,MeanRandom=rMeanWU,summary(m6)$p.table[,1],summary(m6)$s.table[,1],summary(m6)$p.table[,4],summary(m6)$s.table[,4],summary(m6)$dev.expl*100,NbDyads=nrow(dataset1)/2,Type="WU")
      res<-data.frame(rbind(jd1,jd2,jd3,jd4,jd5,jd6))
      names(res)[10:17]<-c("p.(Intercept)","p.ObservedPair","p.SexM","p.TotReads","p.s.AgeMonth","p.s.Ind1","p.s.Ind2","dev.expl") 
      return(res)
    }  
    
    #Models On all pairs
    mod1<-mod(all_pair) 
    all<-rbind(all,mod1)
    
    #Models On young infants only
    mod2<-mod(young_pair)
    young<-rbind(young,mod2)
    
    #Models On old infants only
    mod3<-mod(old_pair)
    old<-rbind(old,mod3)
    
    #Models On old infants only and subsampled
    mod4<-mod(old_pair_sub)
    old_sub1<-rbind(old_sub1,mod4)
    
  }
  
  list <- list()
  list[[1]] <- data.frame(all)
  list[[2]] <- data.frame(young)
  list[[3]] <- data.frame(old)
  list[[4]] <- data.frame(old_sub1)
  return(list)
}


#Generate table results
table_results<-function(res,n_simu){
  res$MeanObs<-as.numeric(as.character(res$MeanObs)) 
  res$MeanRandom<-as.numeric(as.character(res$MeanRandom)) 
  res$ObservedPair1<-as.numeric(as.character(res$ObservedPair1))
  
  nb<-droplevels(subset(res,res$Type %in% "NumberSharedASV"))
  prop.inf<-droplevels(subset(res,res$Type %in% "PropSharedASV.inf"))
  prop.all<-droplevels(subset(res,res$Type %in% "PropSharedASV.all"))
  bc<-droplevels(subset(res,res$Type %in% "BC"))
  uu<-droplevels(subset(res,res$Type %in% "UU"))
  wu<-droplevels(subset(res,res$Type %in% "WU"))
  
  OBS<-c(mean(nb$MeanObs),
         mean(prop.inf$MeanObs),
         mean(prop.all$MeanObs),
         mean(bc$MeanObs),
         mean(uu$MeanObs),
         mean(wu$MeanObs))
  
  SIMU<-c(mean(nb$MeanRandom),
          mean(prop.inf$MeanRandom),
          mean(prop.all$MeanRandom),
          mean(bc$MeanRandom),
          mean(uu$MeanRandom),
          mean(wu$MeanRandom))
  
  quantile<-rbind(quantile(nb$ObservedPair1,probs = c(0.025,0.975)),
                  quantile(prop.inf$ObservedPair1,probs = c(0.025,0.975)),
                  quantile(prop.all$ObservedPair1,probs = c(0.025,0.975)),
                  quantile(bc$ObservedPair1,probs = c(0.025,0.975)),
                  quantile(uu$ObservedPair1,probs = c(0.025,0.975)),
                  quantile(wu$ObservedPair1,probs = c(0.025,0.975)))
  
  pval<-c( 
    length(nb$ObservedPair1[nb$ObservedPair1<0])/n_simu,
    length(prop.inf$ObservedPair1[prop.inf$ObservedPair1<0])/n_simu,
    length(prop.all$ObservedPair1[prop.all$ObservedPair1<0])/n_simu,
    length(bc$ObservedPair1[bc$ObservedPair1>0])/n_simu,
    length(uu$ObservedPair1[uu$ObservedPair1>0])/n_simu,
    length(wu$ObservedPair1[wu$ObservedPair1>0])/n_simu)
  
  tab<-cbind(OBS,SIMU,quantile,pval)
  rownames(tab)<-c("nb","prop.inf","prop.all","bc","uu","wu")
  return(tab)
}

n_simu<-5000
res<-loop(pool3,n_simu) #takes a long long time (it is pre-calculated if needed - see below)
write.csv(res[[1]],"new datasets/mother-offspring pairs/CSS_60days_all pairs.txt")
write.csv(res[[2]],"new datasets/mother-offspring pairs/CSS_60days_young pairs.txt")
write.csv(res[[3]],"new datasets/mother-offspring pairs/CSS_60days_old pairs.txt")
write.csv(res[[4]],"new datasets/mother-offspring pairs/CSS_60days_old pairs subsampled.txt")

table_results(res[[1]],n_simu) #all pairs
table_results(res[[2]],n_simu) #young pairs
table_results(res[[3]],n_simu) #old pairs
table_results(res[[4]],n_simu) #old pairs subsampled

df1<-res[[1]]
df2<-res[[2]]
df3<-res[[3]]
df4<-res[[4]]

# PLOT #
setwd("~/github/")
df1<-read.csv("datasets/mother-offspring pairs/CSS_60days_all pairs.txt",sep=",",header=TRUE)
df2<-read.csv("datasets/mother-offspring pairs/CSS_60days_young pairs.txt",sep=",",header=TRUE)
df3<-read.csv("datasets/mother-offspring pairs/CSS_60days_old pairs.txt",sep=",",header=TRUE)
df4<-read.csv("datasets/mother-offspring pairs/CSS_60days_old pairs subsampled.txt",sep=",",header=TRUE)

table_results(df1,5000) #all pairs
table_results(df2,5000) #young pairs
table_results(df3,5000) #old pairs
table_results(df4,5000) #old pairs subsampled


plot_results<-function(df){
  df$MeanObs<-as.numeric(as.character(df$MeanObs)) 
  df$MeanRandom<-as.numeric(as.character(df$MeanRandom)) 
  prop.inf<-droplevels(subset(df,df$Type %in% "PropSharedASV.inf"))
  prop.all<-droplevels(subset(df,df$Type %in% "PropSharedASV.all"))
  bc<-droplevels(subset(df,df$Type %in% "BC"))
  uu<-droplevels(subset(df,df$Type %in% "UU")) 
  wu<-droplevels(subset(df,df$Type %in% "WU"))
  
  hist(prop.inf$MeanRandom,xlim=c(min(prop.inf$MeanObs,prop.inf$MeanRandom)-1,max(prop.inf$MeanObs,prop.inf$MeanRandom)+1),nclass=100,main="Proportion of shared ASVs (divided inf reads)",xlab="")
  abline(v=mean(prop.inf$MeanObs),col="coral",lwd=4)
  
  hist(prop.all$MeanRandom,xlim=c(min(prop.all$MeanObs,prop.all$MeanRandom)-0.2,max(prop.all$MeanObs,prop.all$MeanRandom)+0.2),nclass=100,main="Proportion of shared ASVs (divided AF+inf reads)",xlab="")
  abline(v=mean(prop.all$MeanObs),col="coral",lwd=4)
  
  hist(bc$MeanRandom,xlim=c(min(bc$MeanObs,bc$MeanRandom)-0.005,max(bc$MeanObs,bc$MeanRandom)+0.005),nclass=100,main="Bray-Curtis",xlab="")
  abline(v=mean(bc$MeanObs),col="coral",lwd=4)
  
  hist(uu$MeanRandom,xlim=c(min(uu$MeanObs,uu$MeanRandom)-0.005,max(uu$MeanObs,uu$MeanRandom)+0.005),nclass=100,font.lab=2,cex.lab=1.2,main="Unweighted Unifrac distance",xlab="")
  abline(v=mean(uu$MeanObs),col="coral",lwd=4)
  
  hist(wu$MeanRandom,nclass=100,xlim=c(min(wu$MeanObs,wu$MeanRandom)-0.0015,max(wu$MeanObs,wu$MeanRandom)+0.001),main="Weighted Unifrac distance",xlab="")
  abline(v=mean(wu$MeanObs),col="coral",lwd=4)
  
}

par(mfcol=c(5,4))
plot_results(df1)
plot_results(df2)
plot_results(df3)
plot_results(df4)
