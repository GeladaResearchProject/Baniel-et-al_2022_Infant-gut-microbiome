rm(list=ls())
'%!in%'<-function(x,y)!('%in%'(x,y))
library(ggplot2)
library(mgcv)  
library(mgcViz)
library(lubridate)
theme_set(theme_classic(base_size=15))

setwd("~/github/")

###### Dataset of compositional similarity (or dissimilarity) between all possible pairs of infant-adult female samples ######
diss<-readRDS("datasets/mother-offspring pairs/CompositionalSimilarity_ALL PAIRS_css_ASV 500 reads.txt")
head(diss) 
dim(diss)#325.488
summary(diss)
diss$Match<-paste(diss$SampleID1,diss$SampleID2,sep="")
diss$Ind1<-as.factor(diss$Ind1)
diss$Ind2<-as.factor(diss$Ind2)
table(diss$SampleType1) #SampleID1 = always an infant sample
table(diss$SampleType2) #SampleID2 = always an adult female sample
diss$TotReads<-diss$NumberReads1+diss$NumberReads2
diss$logTotReads<-log10(diss$TotReads)
diss$ReproState2<-as.character(diss$ReproState2)
diss$ReproState2[diss$ReproState2 %in% "Non-lactating"]<-"aNon-lactating"
diss$ReproState2<-as.factor(diss$ReproState2)

#Add mom ID
pair<-readRDS("datasets/metadata.txt")
head(pair)
dim(pair)
diss1<-merge(diss,pair[c("SampleID","Mother")],by.x="SampleID1",by.y="SampleID",all.x=TRUE,all.y=FALSE)
diss1$MotherSample<-ifelse(as.character(diss1$Ind2)==as.character(diss1$Mother),"1","0")

#Only pairs within 20 days of each other
pool1<-droplevels(subset(diss1,diss1$DiffDay<=20))
dim(pool1)#17.196

#Only pairs from same unit 
pool2<-droplevels(subset(pool1,pool1$SameUnit %in% 1))
dim(pool2)#2319
table(pool2$MotherSample)#558 mom-offspring possible pairs 


#### RESAMPLING ####
loop<-function(pool3,n_simu){
  
  af<-droplevels(subset(pool3,pool3$MotherSample %in% 0))
  mom<-droplevels(subset(pool3,pool3$MotherSample %in% 1))
  good<-intersect(unique(mom$SampleID1),unique(af$SampleID1))
  af<-droplevels(subset(af,af$SampleID1 %in% good))
  mom<-droplevels(subset(mom,mom$SampleID1 %in% good))
  
  
  all<-NULL
  young<-NULL
  old<-NULL
  old_sub1<-NULL
  
  for (j in 1:n_simu) {
    print(j)
    
    ### Generate random pool of females ###
    obs<-NULL
    rand<-NULL
    obs<-sapply(unique(mom$SampleID1),function(x) sample(mom$SampleID2[mom$SampleID1==x],1,replace=FALSE))
    obs1<-data.frame(cbind(SampleID1=as.character(unique(mom$SampleID1)),SampleID2=as.character(obs),MotherSample=rep(1,length(obs))))
    rand<-sapply(unique(af$SampleID1),function(x) sample(af$SampleID2[af$SampleID1==x],1,replace=FALSE))
    rand1<-data.frame(cbind(SampleID1=as.character(unique(af$SampleID1)),SampleID2=as.character(rand),MotherSample=rep(0,length(rand))))
    jd<-rbind(obs1,rand1)
    jd$Match<-paste(jd$SampleID1,jd$SampleID2,sep="")
    dataset<-droplevels(subset(diss1,diss1$Match %in% jd$Match))
    
    all_pair<-dataset
    young_pair<-droplevels(subset(dataset,dataset$AgeMonth1<=12))
    old_pair<-droplevels(subset(dataset,dataset$AgeMonth1>=18))
    
    #Subsample in the old samples to have only same than in young sample
    sub<-sample(unique(old_pair$SampleID1),size=length(unique(young_pair$SampleID1)),replace=FALSE)
    old_pair_sub<-droplevels(subset(old_pair,old_pair$SampleID1 %in% sub))
    
    mod<-function(dataset1){
      res<-NULL
      observed<-droplevels(subset(dataset1,dataset1$MotherSample %in% 1))
      oMeanNumberSharedASV<-mean(observed$NumberSharedASV)
      oMeanPropSharedASV.inf<-mean(observed$PropSharedASV.inf)
      oMeanPropSharedASV.all<-mean(observed$PropSharedASV.all)
      oMeanBC<-mean(observed$bcval)
      oMeanUU<-mean(observed$uval)
      oMeanWU<-mean(observed$wval)
      
      random<-droplevels(subset(dataset1,dataset1$MotherSample %in% 0))
      rMeanNumberSharedASV <-mean(random$NumberSharedASV)
      rMeanPropSharedASV.inf<-mean(random$PropSharedASV.inf)
      rMeanPropSharedASV.all<-mean(random$PropSharedASV.all)
      rMeanBC<-mean(random$bcval)
      rMeanUU<-mean(random$uval)
      rMeanWU<-mean(random$wval)
      
      dataset1$ObservedPair<-factor(dataset1$MotherSample,levels = c("0","1"))
      m1<-gam(NumberSharedASV ~ ObservedPair + s(AgeMonth1) + ReproState2 + DiffDay + Sex1 + logTotReads + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset1,gamma=1.4,method="REML") 
      m2<-gam(PropSharedASV.inf ~ ObservedPair + s(AgeMonth1)  + ReproState2 + DiffDay + Sex1 + logTotReads + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset1,gamma=1.4,method="REML") 
      m3<-gam(PropSharedASV.all ~ ObservedPair + s(AgeMonth1)  + ReproState2 + DiffDay + Sex1 + logTotReads + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset1,gamma=1.4,method="REML") 
      m4<-gam(bcval ~ ObservedPair + s(AgeMonth1)  + ReproState2 + DiffDay + Sex1  + logTotReads + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset1,gamma=1.4,method="REML") 
      m5<-gam(uval ~ ObservedPair + s(AgeMonth1)  + ReproState2 + DiffDay + Sex1 + logTotReads + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset1,gamma=1.4,method="REML") 
      m6<-gam(wval ~ ObservedPair+ s(AgeMonth1)  + ReproState2 + DiffDay + Sex1 + logTotReads + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset1,gamma=1.4,method="REML") 
      
      jd1<-c(MeanObs=oMeanNumberSharedASV,MeanRandom=rMeanNumberSharedASV,summary(m1)$p.table[,1],summary(m1)$s.table[,1],summary(m1)$p.table[,4],summary(m1)$s.table[,4],summary(m1)$dev.expl*100,NbDyads=nrow(dataset1)/2,Type="NumberSharedASV")
      jd2<-c(MeanObs=oMeanPropSharedASV.inf,MeanRandom=rMeanPropSharedASV.inf,summary(m2)$p.table[,1],summary(m2)$s.table[,1],summary(m2)$p.table[,4],summary(m2)$s.table[,4],summary(m2)$dev.expl*100,NbDyads=nrow(dataset1)/2,Type="PropSharedASV.inf")
      jd3<-c(MeanObs=oMeanPropSharedASV.all,MeanRandom=rMeanPropSharedASV.all,summary(m3)$p.table[,1],summary(m3)$s.table[,1],summary(m3)$p.table[,4],summary(m3)$s.table[,4],summary(m3)$dev.expl*100,NbDyads=nrow(dataset1)/2,Type="PropSharedASV.all")
      jd4<-c(MeanObs=oMeanBC,MeanRandom=rMeanBC,summary(m4)$p.table[,1],summary(m4)$s.table[,1],summary(m4)$p.table[,4],summary(m4)$s.table[,4],summary(m4)$dev.expl*100,NbDyads=nrow(dataset1)/2,Type="BC")
      jd5<-c(MeanObs=oMeanUU,MeanRandom=rMeanUU,summary(m5)$p.table[,1],summary(m5)$s.table[,1],summary(m5)$p.table[,4],summary(m5)$s.table[,4],summary(m5)$dev.expl*100,NbDyads=nrow(dataset1)/2,Type="UU")
      jd6<-c(MeanObs=oMeanWU,MeanRandom=rMeanWU,summary(m6)$p.table[,1],summary(m6)$s.table[,1],summary(m6)$p.table[,4],summary(m6)$s.table[,4],summary(m6)$dev.expl*100,NbDyads=nrow(dataset1)/2,Type="WU")
      res<-data.frame(rbind(jd1,jd2,jd3,jd4,jd5,jd6))
      names(res)[12:21]<-c("p.(Intercept)","p.ObservedPair","p.ReproStateL","p. DiffDay","p.SexM","p.logTotReads","p.s.AgeMonth","p.s.Ind1","p.s.Ind2","dev.expl") 
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
  bc<-droplevels(subset(res,res$Type %in% "BC"))
  uu<-droplevels(subset(res,res$Type %in% "UU"))
  wu<-droplevels(subset(res,res$Type %in% "WU"))
  
  OBS<-c(mean(nb$MeanObs),
         mean(bc$MeanObs),
         mean(uu$MeanObs),
         mean(wu$MeanObs))
  
  SIMU<-c(mean(nb$MeanRandom),
          mean(bc$MeanRandom),
          mean(uu$MeanRandom),
          mean(wu$MeanRandom))
  
  estimate<-rbind(mean(nb$ObservedPair1),
                  mean(bc$ObservedPair1),
                  mean(uu$ObservedPair1),
                  mean(wu$ObservedPair1))
  
  quantile<-rbind(quantile(nb$ObservedPair1,probs = c(0.025,0.975)),
                  quantile(bc$ObservedPair1,probs = c(0.025,0.975)),
                  quantile(uu$ObservedPair1,probs = c(0.025,0.975)),
                  quantile(wu$ObservedPair1,probs = c(0.025,0.975)))
  
  pval<-c( 
    length(nb$ObservedPair1[nb$ObservedPair1<0])/n_simu,
    length(bc$ObservedPair1[bc$ObservedPair1>0])/n_simu,
    length(uu$ObservedPair1[uu$ObservedPair1>0])/n_simu,
    length(wu$ObservedPair1[wu$ObservedPair1>0])/n_simu)
  
  tab<-cbind(OBS,SIMU,estimate,quantile,pval)
  rownames(tab)<-c("nb","bc","uu","wu")
  return(tab)
}


n_simu<-5000
res<-loop(pool2,n_simu) #takes a long long time (it is pre-calculated if needed - see below)
write.csv(res[[1]],"datasets/mother-offspring pairs/CSS_20days_all pairs.txt")
write.csv(res[[2]],"datasets/mother-offspring pairs/CSS_20days_young pairs.txt")
write.csv(res[[3]],"datasets/mother-offspring pairs/CSS_20days_old pairs.txt")
write.csv(res[[4]],"datasets/mother-offspring pairs/CSS_20days_old pairs subsampled.txt")


# RESULTS & PLOT #
setwd("~/github/datasets/mother-offspring pairs/")
df1<-read.csv("CSS_20days_all pairs.txt",sep=",",header=TRUE)
df2<-read.csv("CSS_20days_young pairs.txt",sep=",",header=TRUE)
df3<-read.csv("CSS_20days_old pairs.txt",sep=",",header=TRUE)
df4<-read.csv("CSS_20days_old pairs subsampled.txt",sep=",",header=TRUE)

table_results(df1,5000) #all pairs
table_results(df2,5000) #young pairs
table_results(df3,5000) #old pairs
table_results(df4,5000) #old pairs subsampled


plot_results<-function(df,main){
  df$MeanObs<-as.numeric(as.character(df$MeanObs)) 
  df$MeanRandom<-as.numeric(as.character(df$MeanRandom)) 
  nb<-droplevels(subset(df,df$Type %in% "NumberSharedASV"))
  bc<-droplevels(subset(df,df$Type %in% "BC"))
  uu<-droplevels(subset(df,df$Type %in% "UU")) 
  wu<-droplevels(subset(df,df$Type %in% "WU"))
  
  hist(nb$MeanRandom,xlim=c(min(nb$MeanObs,nb$MeanRandom)-10,max(nb$MeanObs,nb$MeanRandom)+10),nclass=100,cex.lab=1.2,xlab="# shared ASVs",main=main)
  abline(v=mean(nb$MeanObs),col="coral",lwd=4)
  
  hist(bc$MeanRandom,xlim=c(min(bc$MeanObs,bc$MeanRandom)-0.005,max(bc$MeanObs,bc$MeanRandom)+0.005),nclass=100,cex.lab=1.2,xlab="Bray-Curtis",main=main)
  abline(v=mean(bc$MeanObs),col="coral",lwd=4)
  
  hist(uu$MeanRandom,xlim=c(min(uu$MeanObs,uu$MeanRandom)-0.005,max(uu$MeanObs,uu$MeanRandom)+0.005),nclass=100,cex.lab=1.2,xlab="Unweighted Unifrac",main=main)
  abline(v=mean(uu$MeanObs),col="coral",lwd=4)
  
  hist(wu$MeanRandom,nclass=100,xlim=c(min(wu$MeanObs,wu$MeanRandom)-0.0015,max(wu$MeanObs,wu$MeanRandom)+0.001),cex.lab=1.2,xlab="Weighted Unifrac",main=main)
  abline(v=mean(wu$MeanObs),col="coral",lwd=4)
  
}

par(mfcol=c(4,4))
plot_results(df1,"All immatures")
plot_results(df2,"Young infants")
plot_results(df3,"Old juveniles")
plot_results(df4,"Old juveniles (subset)")






