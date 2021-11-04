rm(list=ls())
'%!in%'<-function(x,y)!('%in%'(x,y))
library(ggplot2)
library(mgcv)  
library(mgcViz)
theme_set (theme_classic(base_size=15))

setwd("~/github/")

###### Dataset of compositional similarity (or dissimilarity) between all possible pairs of infant-adult female samples ######
diss<-readRDS("datasets/mother-offspring pairs/CompositionalSimilarity_ALL PAIRS.txt")
head(diss) 
dim(diss)
diss$Ind1<-as.factor(diss$Ind1)
diss$Ind2<-as.factor(diss$Ind2)
#SampleID1 = always an infant sample, SampleID2 = always an adult female sample
#PropSharedASV = number of shared ASV between the infant and female samples divided by number of ASVs in the infant sample 

### Observed similarity (between 398 actual mother-offspring matched pairs) ###  
pair<-readRDS("datasets/metadata.txt")
pair<-droplevels(subset(pair,pair$MatchedMomSample %!in% NA))
pair$Match<-paste(pair$SampleID,pair$MatchedMomSample,sep="")
head(pair)
dim(pair) #398 mother-infant pairs (collected0 or 1 day apart) 

observed<-merge(pair,diss[c("Match","SameSeason","SameUnit","bcval","uval","wval","NumberASV1","NumberASV2","NumberSharedASV","PropSharedASV")],by="Match",all.x=TRUE,all.y=FALSE)
head(observed)
summary(observed) 

#AGE 
ggplot(observed, aes(x=AgeMonth,y=NumberSharedASV,col=Weaned)) + geom_point() + xlab("Age (month)") + ylab("Number of shared ASV between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))
ggplot(observed, aes(x=AgeMonth,y=PropSharedASV,col=Weaned)) + geom_point() + xlab("Age (month)") + ylab("Proportion of shared ASV between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))
ggplot(observed, aes(x=AgeMonth,y=bcval,col=Weaned)) + geom_point() + xlab("Age (months)") + ylab("Bray-Curtis dissimilarity between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))
ggplot(observed, aes(x=AgeMonth,y=uval,col=Weaned)) + geom_point() + xlab("Age (months)") + ylab("Unweigthed Unifrac dissimilarity between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))
ggplot(observed, aes(x=AgeMonth,y=wval,col=Weaned)) + geom_point() + xlab("Age (months)") + ylab("Weigthed Unifrac dissimilarity between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))

#PARTIY
ggplot(observed, aes(x=AgeMonth,y=NumberSharedASV,col=Parity)) + geom_point() + xlab("Age (month)") + ylab("Number of shared ASV between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)
ggplot(observed, aes(x=AgeMonth,y=PropSharedASV,col=Parity)) + geom_point() + xlab("Age (month)") + ylab("Proportion of shared ASV between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)
ggplot(observed, aes(x=AgeMonth,y=bcval,col=Parity)) + geom_point() + xlab("Age (months)") + ylab("Bray-Curtis dissimilarity between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)
ggplot(observed, aes(x=AgeMonth,y=uval,col=Parity)) + geom_point() + xlab("Age (months)") + ylab("Unweigthed Unifrac dissimilarity between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)
ggplot(observed, aes(x=AgeMonth,y=wval,col=Parity)) + geom_point() + xlab("Age (months)") + ylab("Weigthed Unifrac dissimilarity between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)


#GAMMs - ALL INFANTS
head(observed)
dim(observed) #398 pairs

m1<-gam(NumberSharedASV ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30, data=observed,gamma=1.4,method="REML")
summary(m1)

m2<-gam(uval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30, data=observed,gamma=1.4,method="REML")
summary(m2)

m3<-gam(wval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30, data=observed,gamma=1.4,method="REML")
summary(m3)


#GAMMs - YOUNG INFANTS
young<-droplevels(subset(observed,observed$AgeMonth<=12))
dim(young) #136 pairs 

m1<-gam(NumberSharedASV ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re")  + Sex + Parity + MomRank + Rain30 + MinT30, data=young,gamma=1.4,method="REML")
summary(m1)
ggplot(young, aes(x=AgeMonth,y=NumberSharedASV,col=Parity)) + geom_point() + xlab("Age (month)") + ylab("Number of shared ASV between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)

m2<-gam(uval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30, data=young,gamma=1.4,method="REML")
summary(m2)
ggplot(young, aes(x=AgeMonth,y=uval,col=Parity)) + geom_point() + xlab("Age (month)") + ylab("Number of shared ASV between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)

m3<-gam(wval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30, data=young,gamma=1.4,method="REML")
summary(m3)
ggplot(young, aes(x=AgeMonth,y=wval,col=Parity)) + geom_point() + xlab("Age (month)") + ylab("Number of shared ASV between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)
ggplot(young, aes(x=AgeMonth,y=wval,col=MomRank_cat)) + geom_point() + xlab("Age (month)") + ylab("Number of shared ASV between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)


#GAMMs - OLD INFANTS
old<-droplevels(subset(observed,observed$AgeMonth>=18))
dim(young) #136 pairs

m1<-gam(NumberSharedASV ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re")  + Sex + Parity + MomRank + Rain30 + MinT30, data=old,gamma=1.4,method="REML")
summary(m1)

m2<-gam(uval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30, data=old,gamma=1.4,method="REML")
summary(m2)

m3<-gam(wval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30, data=old,gamma=1.4,method="REML")
summary(m3)


### QUADRATIC PLATEAU MODELS ###
library(easynls)
data=observed[c("AgeMonth","NumberSharedASV")]
nlsplot (data,model=4,xlab="Age",ylab="NumberSharedASV")
nlsfit(data,model=4)
dta.nls <- nls(NumberSharedASV ~ (a + b * AgeMonth + c * I(AgeMonth^2)) * (AgeMonth <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (AgeMonth > -0.5 * b/c), observed,
               start=list(a=35.55076462,b=as.numeric(51.45726147),c=as.numeric(-1.77704424)))
observed$line <- predict(dta.nls)
z<-ggplot(observed,aes(x=AgeMonth, y=NumberSharedASV,colour=Weaned))+ geom_point(size=2)+ 
  scale_color_manual(values=c("#440154FF", "#5DC863FF"))  +
  xlab("Age (months)") + ylab("Number of shared ASVs with mother")+ 
  theme(legend.position="none",axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  geom_line(aes(y =line),size=1.8,color='coral')+
  annotate("segment", x = 14.47832877, 
           xend = 14.47832877, y = 0, 
           yend = 400,colour = "coral",size=1,linetype="dashed")+
  annotate(geom="text", x=25, y=100, label="14.48 months",size=4,color="coral")
z1<-z+ labs(color='Immature weaned') +  theme(legend.title = element_text(size=15,face="bold"),legend.text =element_text(size=15))
z1


data=observed[c("AgeMonth","uval")]
nlsplot (data,model=4,xlab="Age",ylab="Observed")
nlsfit(data,model=4)
dta.nls <- nls(uval ~ (a + b * AgeMonth + c * I(AgeMonth^2)) * (AgeMonth <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (AgeMonth > -0.5 * b/c), observed,
               start=list(a=0.89468828,b=as.numeric(-0.05197348),c=as.numeric( 0.00167597)))
observed$line <- predict(dta.nls)
w<-ggplot(observed,aes(x=AgeMonth, y=uval,colour=Weaned))+ geom_point(size=2)+ 
  scale_color_manual(values=c("#440154FF", "#5DC863FF"))  +
  xlab("Age (months)") + ylab("Unweigthed Unifrac dissimilarity")+ 
  theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"))+
  geom_line(aes(y =line),size=1.8,color='coral')+
  annotate("segment", x = 15.50545445, 
           xend = 15.50545445, y = 0.49, 
           yend = 0.94,colour = "coral",size=1,linetype="dashed")+
  annotate(geom="text", x=25, y=0.95, label="15.51 months",size=4,color="coral")
w1<-w+ labs(color='Immature weaned') +  theme(legend.title = element_text(size=15,face="bold"),legend.text =element_text(size=15))
w1


data=observed[c("AgeMonth","wval")]
nlsplot (data,model=4,xlab="Age",ylab="Observed")
nlsfit(data,model=4)
dta.nls <- nls(wval ~ (a + b * AgeMonth + c * I(AgeMonth^2)) * (AgeMonth <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (AgeMonth > -0.5 * b/c), observed,
               start=list(a=1.010416e-01,b=as.numeric( -5.200850e-03),c=as.numeric(1.761300e-04)))
observed$line <- predict(dta.nls)

q<-ggplot(observed,aes(x=AgeMonth, y=wval,colour=Weaned))+ geom_point(size=2)+ 
  scale_color_manual(values=c("#440154FF", "#5DC863FF"))  +
  xlab("Age (months)") + ylab("Weigthed Unifrac dissimilarity")+ 
  theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"))+
  geom_line(aes(y =line),size=1.8,color='coral')+
  annotate("segment", x = 1.476403e+01, 
           xend = 1.476403e+01, y = 0.065, 
           yend = 0.18,colour = "coral",size=1,linetype="dashed")+
  annotate(geom="text", x=25, y=0.19, label="14.76 months",size=4,color="coral")
q1<-q+ labs(color='Immature weaned') +  theme(legend.title = element_text(size=15,face="bold"),legend.text =element_text(size=15))
q1






################################# RESAMPLING APPROACH ################################# 
#Test if mother-offspring pairs are more similar in composition than expected by chance (i.e, with a random female of the population) 

#Add info of the actual mother matched sample
head(diss)
diss1<-merge(diss,pair[c("SampleID","Mother","MatchedMomSample")],by.x="SampleID1",by.y="SampleID",all.x=TRUE,all.y=FALSE)
diss1<-diss1[c(2,1,3,4:28)]
head(diss1)
names(diss1)
summary(diss)

diss1$ObservedPair<-ifelse(as.character(diss1$SampleID2)==as.character(diss1$MatchedMomSample),"1","0")
table(diss1$ObservedPair)
diss1$ObservedPair<-as.factor(diss1$ObservedPair)

diss1$MotherSample<-ifelse(as.character(diss1$Ind2)==as.character(diss1$Mother),"1","0")
table(diss1$MotherSample)

#Observed dataset 
observed<-droplevels(subset(diss1,diss1$ObservedPair %in% "1")) 
dim(observed)

observed_young<-droplevels(subset(observed,observed$AgeMonth1<=12))
dim(observed_young)

observed_old<-droplevels(subset(observed,observed$AgeMonth1>=18))
dim(observed_old)

#Create pool of random females 
pool<-droplevels(subset(diss1,diss1$MotherSample %!in% "1")) #remove all matched mother samples from the pool of random females
table(pool$ObservedPair)#no true ObservedPair anymore in the random pool
dim(pool)#240429
summary(pool)



#### RESAMPLING ####
loop<-function(pool_rand,n_simu){
  
all<-NULL
young<-NULL
old<-NULL

for (j in 1:n_simu) {
  print(j)

### Generate random pool of females ###
rand<-NULL
rand<-sapply(unique(pool_rand$SampleID1),function(x) sample(pool_rand$SampleID2[pool_rand$SampleID1==x],1,replace=FALSE))
rand1<-data.frame(cbind(SampleID1=as.character(unique(pool_rand$SampleID1)),SampleID2=as.character(rand)))
rand1$Match<-paste(rand1$SampleID1,rand1$SampleID2,sep="")
rand2<-merge(rand1,diss1[,-c(2:3)],by="Match",all.x=TRUE,all.y=FALSE)
dataset<-rbind(observed,rand2)

mod<-function(dataset){
  res<-NULL
  random<-droplevels(subset(dataset,dataset$ObservedPair %in% 0))
  MeanNumberSharedASV<-mean(random$NumberSharedASV)
  MeanPropSharedASV<-mean(random$PropSharedASV)
  MeanBC<-mean(random$bcval)
  MeanUU<-mean(random$uval)
  MeanWU<-mean(random$wval)
  dataset$ObservedPair<-factor(dataset$ObservedPair,levels = c("0","1"))
  m1<-gam(NumberSharedASV ~ ObservedPair + s(AgeMonth1) + Sex1 + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset,gamma=1.4,method="REML") 
  m2<-gam(PropSharedASV ~ ObservedPair + s(AgeMonth1) + Sex1 + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset,gamma=1.4,method="REML") 
  m3<-gam(bcval ~ ObservedPair + s(AgeMonth1) + Sex1  + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset,gamma=1.4,method="REML") 
  m4<-gam(uval ~ ObservedPair + s(AgeMonth1) + Sex1 + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset,gamma=1.4,method="REML") 
  m5<-gam(wval ~ ObservedPair+ s(AgeMonth1) + Sex1 + s(Ind1,bs="re")+ s(Ind2,bs="re"), data=dataset,gamma=1.4,method="REML") 
  
  jd1<-c(MeanRandom=MeanNumberSharedASV,summary(m1)$p.table[,1],summary(m1)$s.table[,1],summary(m1)$p.table[,4],summary(m1)$s.table[,4],summary(m1)$dev.expl*100,NbDyads=nrow(dataset)/2,Type="NumberSharedASV")
  jd2<-c(MeanRandom=MeanPropSharedASV,summary(m2)$p.table[,1],summary(m2)$s.table[,1],summary(m2)$p.table[,4],summary(m2)$s.table[,4],summary(m2)$dev.expl*100,NbDyads=nrow(dataset)/2,Type="PropSharedASV")
  jd3<-c(MeanRandom=MeanBC,summary(m3)$p.table[,1],summary(m3)$s.table[,1],summary(m3)$p.table[,4],summary(m3)$s.table[,4],summary(m3)$dev.expl*100,NbDyads=nrow(dataset)/2,Type="BC")
  jd4<-c(MeanRandom=MeanUU,summary(m4)$p.table[,1],summary(m4)$s.table[,1],summary(m4)$p.table[,4],summary(m4)$s.table[,4],summary(m4)$dev.expl*100,NbDyads=nrow(dataset)/2,Type="UU")
  jd5<-c(MeanRandom=MeanWU,summary(m5)$p.table[,1],summary(m5)$s.table[,1],summary(m5)$p.table[,4],summary(m5)$s.table[,4],summary(m5)$dev.expl*100,NbDyads=nrow(dataset)/2,Type="WU")
  res<-data.frame(rbind(jd1,jd2,jd3,jd4,jd5))
  names(res)[8:14]<-c("p.(Intercept)","p.ObservedPair","p.SexM","p.s.AgeMonth","p.s.Ind1","p.s.Ind2","dev.expl") 
  return(res)
}  

#Models On all pairs
mod1<-mod(dataset) 
all<-rbind(all,mod1)
  
#Models On young infants only
young_pair<-droplevels(subset(dataset,dataset$AgeMonth1<=12))
mod2<-mod(young_pair)
young<-rbind(young,mod2)
  
#Models On old infants only
old_pair<-droplevels(subset(dataset,dataset$AgeMonth1>18))
mod3<-mod(old_pair)
old<-rbind(old,mod3)
}

list <- list()
list[[1]] <- data.frame(all)
list[[2]] <- data.frame(young)
list[[3]] <- data.frame(old)

return(list)
}


#Generate table results
table_results<-function(res,obs,n_simu){
  res$MeanRandom<-as.numeric(as.character(res$MeanRandom)) 
  res$ObservedPair1<-as.numeric(as.character(res$ObservedPair1))
  
  nb<-droplevels(subset(res,res$Type %in% "NumberSharedASV"))
  prop<-droplevels(subset(res,res$Type %in% "PropSharedASV"))
  bc<-droplevels(subset(res,res$Type %in% "BC"))
  uu<-droplevels(subset(res,res$Type %in% "UU"))
  wu<-droplevels(subset(res,res$Type %in% "WU"))
  
  OBS<-c(mean(obs$NumberSharedASV),
         mean(obs$PropSharedASV),
         mean(obs$bcval),
         mean(obs$uval),
         mean(obs$wval))
  
  SIMU<-c(mean(nb$MeanRandom),
          mean(prop$MeanRandom),
          mean(bc$MeanRandom),
          mean(uu$MeanRandom),
          mean(wu$MeanRandom))
  
  quantile<-rbind(quantile(nb$ObservedPair1,probs = c(0.025,0.975)),
                  quantile(prop$ObservedPair1,probs = c(0.025,0.975)),
                  quantile(bc$ObservedPair1,probs = c(0.025,0.975)),
                  quantile(uu$ObservedPair1,probs = c(0.025,0.975)),
                  quantile(wu$ObservedPair1,probs = c(0.025,0.975)))
  
  pval<-c(
    length(nb$ObservedPair1[nb$ObservedPair1<0])/n_simu,
    length(prop$ObservedPair1[prop$ObservedPair1<0])/n_simu,
    length(bc$ObservedPair1[bc$ObservedPair1>0])/n_simu,
    length(uu$ObservedPair1[uu$ObservedPair1>0])/n_simu,
    length(wu$ObservedPair1[wu$ObservedPair1>0])/n_simu)
  
  tab<-cbind(OBS,SIMU,quantile,pval)
  rownames(tab)<-c("nb","prop","bc","uu","wu")
  return(tab)
}


########################## MATCH RANDOM PAIRS BY SEASON ##########################
pool_Season<-droplevels(subset(pool,pool$SameSeason %in% 1))
head(pool_Season)
dim(pool_Season)#83441
summary(as.numeric(table(pool_Season$SampleID1)))#209 possible random match per infant sample on average
n_simu<-1000
#res_Season<-loop(pool_Season,n_simu) #takes a long long time (it is pre-calculated if needed - see below)
#write.csv(res_Season[[1]],"datasets/Match.by.Season_all pairs.txt")
#write.csv(res_Season[[2]],"datasets/Match.by.Season_young pairs.txt")
#write.csv(res_Season[[3]],"datasets/Match.by.Season_old pairs.txt")
#table_results(res_Season[[1]],observed,n_simu) #all pairs
#table_results(res_Season[[2]],observed_young,n_simu) #young pairs
#table_results(res_Season[[3]],observed_old,n_simu) #old pairs


# PLOT #
df<-read.csv("datasets/mother-offspring pairs/Match.by.Season_all pairs.txt",sep=",",header=TRUE); obs=observed
#df<-read.csv("datasets/mother-offspring pairs/Match.by.Season_young pairs.txt",sep=",",header=TRUE); obs=observed_young 
#df<-read.csv("datasets/mother-offspring pairs/Match.by.Season_young pairs.txt",sep=",",header=TRUE); obs=observed_old  

head(df)
df$MeanRandom<-as.numeric(as.character(df$MeanRandom)) 
nb<-droplevels(subset(df,df$Type %in% "NumberSharedASV"))
prop<-droplevels(subset(df,df$Type %in% "PropSharedASV"))
uu<-droplevels(subset(df,df$Type %in% "UU")) #unweighted UniFrac
wu<-droplevels(subset(df,df$Type %in% "WU"))#weighted UniFrac


par(mfrow=c(2,2))

hist(nb$MeanRandom,xlim=c(min(mean(obs$NumberSharedASV),nb$MeanRandom)-1,max(mean(obs$NumberSharedASV),nb$MeanRandom)+1),nclass=100,main="Number of shared ASVs",xlab="")
abline(v=mean(obs$NumberSharedASV),col="coral",lwd=4)

hist(prop$MeanRandom,xlim=c(min(mean(obs$PropSharedASV),prop$MeanRandom)-1,max(mean(obs$PropSharedASV),prop$MeanRandom)+1),nclass=100,main="Proportion of shared ASVs",xlab="")
abline(v=mean(obs$PropSharedASV),col="coral",lwd=4)

hist(uu$MeanRandom,xlim=c(min(mean(obs$uval),uu$MeanRandom)-0.005,max(mean(obs$uval),uu$MeanRandom)+0.005),nclass=100,font.lab=2,cex.lab=1.2,main="Unweighted Unifrac distance",xlab="")
abline(v=mean(obs$uval),col="coral",lwd=4)

hist(wu$MeanRandom,nclass=100,xlim=c(min(mean(obs$wval),wu$MeanRandom)-0.0015,max(mean(obs$wval),wu$MeanRandom)+0.001),main="Weighted Unifrac distance",xlab="")
abline(v=mean(obs$wval),col="coral",lwd=4)



########################## MATCH RANDOM PAIRS BY UNIT ##########################
pool_Unit<-droplevels(subset(pool,pool$SameUnit %in% 1))
head(pool_Unit)
dim(pool_Unit)#22876
summary(as.numeric(table(pool_Unit$SampleID1)))#57 possible random match per infant sample on average
n_simu<-1000
#res_Unit<-loop(pool_Unit,n_simu) ##takes a long long time (it is pre-calculated if needed - see below)
#write.csv(res_Unit[[1]],"datasets/mother-offspring pairs/Match.by.Unit_all pairs.txt")
#write.csv(res_Unit[[2]],"datasets/mother-offspring pairs/Match.by.Unit_young pairs.txt")
#write.csv(res_Unit[[3]],"datasets/mother-offspring pairs/Match.by.Unit_old pairs.txt")
#table_results(res_Unit[[1]],observed,n_simu) #all pairs
#table_results(res_Unit[[2]],observed_young,n_simu) #young pairs
#table_results(res_Unit[[3]],observed_old,n_simu) #old pairs

# PLOT #
df<-read.csv("datasets/mother-offspring pairs/Match.by.Unit_all pairs.txt",sep=",",header=TRUE); obs=observed
#df<-read.csv("datasets/mother-offspring pairs/Match.by.Unit_young pairs.txt",sep=",",header=TRUE); obs=observed_young 
#df<-read.csv("datasets/mother-offspring pairs/Match.by.Unit_young pairs.txt",sep=",",header=TRUE); obs=observed_old  

head(df)
df$MeanRandom<-as.numeric(as.character(df$MeanRandom)) 
nb<-droplevels(subset(df,df$Type %in% "NumberSharedASV"))
prop<-droplevels(subset(df,df$Type %in% "PropSharedASV"))
uu<-droplevels(subset(df,df$Type %in% "UU")) #unweighted UniFrac
wu<-droplevels(subset(df,df$Type %in% "WU"))#weighted UniFrac



par(mfrow=c(2,2))

hist(nb$MeanRandom,xlim=c(min(mean(obs$NumberSharedASV),nb$MeanRandom)-1,max(mean(obs$NumberSharedASV),nb$MeanRandom)+1),nclass=100,main="Number of shared ASVs",xlab="")
abline(v=mean(obs$NumberSharedASV),col="coral",lwd=4)

hist(prop$MeanRandom,xlim=c(min(mean(obs$PropSharedASV),prop$MeanRandom)-1,max(mean(obs$PropSharedASV),prop$MeanRandom)+1),nclass=100,main="Proportion of shared ASVs",xlab="")
abline(v=mean(obs$PropSharedASV),col="coral",lwd=4)

hist(uu$MeanRandom,xlim=c(min(mean(obs$uval),uu$MeanRandom)-0.005,max(mean(obs$uval),uu$MeanRandom)+0.005),nclass=100,font.lab=2,cex.lab=1.2,main="Unweighted Unifrac distance",xlab="")
abline(v=mean(obs$uval),col="coral",lwd=4)

hist(wu$MeanRandom,nclass=100,xlim=c(min(mean(obs$wval),wu$MeanRandom)-0.001,max(mean(obs$wval),wu$MeanRandom)+0.001),main="Weighted Unifrac distance",xlab="")
abline(v=mean(obs$wval),col="coral",lwd=4)

