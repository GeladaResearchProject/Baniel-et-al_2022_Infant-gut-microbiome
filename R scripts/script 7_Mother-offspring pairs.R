rm(list=ls())
'%!in%'<-function(x,y)!('%in%'(x,y))
library(ggplot2)
library(mgcv)  
library(mgcViz)
theme_set (theme_classic(base_size=15))

setwd("~/github/")

###### Dataset of compositional similarity (or dissimilarity) between all possible pairs of infant-adult female samples ######
diss<-readRDS("datasets/mother-offspring pairs/CompositionalSimilarity_ALL PAIRS_css_ASV 500 reads.txt")
head(diss) 
dim(diss)
diss$Ind1<-as.factor(diss$Ind1)
diss$Ind2<-as.factor(diss$Ind2)
table(diss$SampleType1) #SampleID1 = always an infant sample
table(diss$SampleType2) #SampleID2 = always an adult female sample
diss$TotReads<-diss$NumberReads1+diss$NumberReads2

### Observed similarity (between 398 actual mother-offspring matched pairs) ###  
pair<-readRDS("datasets/metadata.txt")
pair<-droplevels(subset(pair,pair$MatchedMomSample %!in% NA))
pair$Match<-paste(pair$SampleID,pair$MatchedMomSample,sep="")
head(pair)
dim(pair) #398 mother-infant pairs (collected0 or 1 day apart) 

observed<-merge(pair,diss[c("Match","SameUnit","ReproState2","NumberReads2","TotReads","DiffRead","DiffDay","bcval","uval","wval","NumberASV1","NumberASV2","NumberSharedASV","PropSharedASV.inf","PropSharedASV.all")],by="Match",all.x=TRUE,all.y=FALSE)
head(observed)
summary(observed) 
colnames(observed)[27:28]<-c("ReproStateMom","NumberReadsMom")
table(observed$SampleType)

#AGE 
ggplot(observed, aes(x=AgeMonth,y=NumberSharedASV,col=Weaned)) + geom_point() + xlab("Age (month)") + ylab("Number of shared ASV between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))
ggplot(observed, aes(x=AgeMonth,y=PropSharedASV.inf,col=Weaned)) + geom_point() + xlab("Age (month)") + ylab("Proportion of shared ASV between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))
ggplot(observed, aes(x=AgeMonth,y=PropSharedASV.all,col=Weaned)) + geom_point() + xlab("Age (month)") + ylab("Proportion of shared ASV between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))
ggplot(observed, aes(x=AgeMonth,y=bcval,col=Weaned)) + geom_point() + xlab("Age (months)") + ylab("Bray-Curtis dissimilarity between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))
ggplot(observed, aes(x=AgeMonth,y=uval,col=Weaned)) + geom_point() + xlab("Age (months)") + ylab("Unweigthed Unifrac dissimilarity between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))
ggplot(observed, aes(x=AgeMonth,y=wval,col=Weaned)) + geom_point() + xlab("Age (months)") + ylab("Weigthed Unifrac dissimilarity between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))

#PARTIY
ggplot(observed, aes(x=AgeMonth,y=NumberSharedASV,col=Parity)) + geom_point() + xlab("Age (month)") + ylab("Number of shared ASV between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)
ggplot(observed, aes(x=AgeMonth,y=PropSharedASV.inf,col=Parity)) + geom_point() + xlab("Age (month)") + ylab("Proportion of shared ASV between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)
ggplot(observed, aes(x=AgeMonth,y=PropSharedASV.all,col=Parity)) + geom_point() + xlab("Age (month)") + ylab("Proportion of shared ASV between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)
ggplot(observed, aes(x=AgeMonth,y=bcval,col=Parity)) + geom_point() + xlab("Age (months)") + ylab("Bray-Curtis dissimilarity between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)
ggplot(observed, aes(x=AgeMonth,y=uval,col=Parity)) + geom_point() + xlab("Age (months)") + ylab("Unweigthed Unifrac dissimilarity between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)
ggplot(observed, aes(x=AgeMonth,y=wval,col=Parity)) + geom_point() + xlab("Age (months)") + ylab("Weigthed Unifrac dissimilarity between mother-offspring pairs")+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) + stat_smooth(method = "gam",alpha=0.2, formula = y ~ s(x), size = 1)


#GAMMs - ALL INFANTS
head(observed)
dim(observed) #398 pairs

m1<-gam(NumberSharedASV ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30 + TotReads, data=observed,gamma=1.4,method="REML")
summary(m1)

m4<-gam(bcval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30 + TotReads, data=observed ,gamma=1.4,method="REML")
summary(m4)

m5<-gam(uval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30 + TotReads, data=observed,gamma=1.4,method="REML")
summary(m5)

m6<-gam(wval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30+ TotReads, data=observed,gamma=1.4,method="REML")
summary(m6)


#GAMMs - YOUNG INFANTS
young<-droplevels(subset(observed,observed$AgeMonth<=12))
dim(young) #136 pairs 

m1<-gam(NumberSharedASV ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30 + TotReads, data=young,gamma=1.4,method="REML")
summary(m1)

m4<-gam(bcval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30 + TotReads, data=young,gamma=1.4,method="REML")
summary(m4)

m5<-gam(uval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30 + TotReads, data=young,gamma=1.4,method="REML")
summary(m5)

m6<-gam(wval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30+ TotReads, data=young,gamma=1.4,method="REML")
summary(m6)


#GAMMs - OLD INFANTS
old<-droplevels(subset(observed,observed$AgeMonth>=18))
dim(young) #136 pairs

m1<-gam(NumberSharedASV ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30 + TotReads, data=old,gamma=1.4,method="REML")
summary(m1)

m4<-gam(bcval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30 + TotReads, data=old,gamma=1.4,method="REML")
summary(m4)

m5<-gam(uval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30 + TotReads, data=old,gamma=1.4,method="REML")
summary(m5)

m6<-gam(wval ~ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re") + Sex + Parity + MomRank + Rain30 + MinT30+ TotReads, data=old,gamma=1.4,method="REML")
summary(m6)



### QUADRATIC PLATEAU MODELS ###
library(easynls)
data=observed[c("AgeMonth","NumberSharedASV")]
nlsplot (data,model=4,xlab="Age",ylab="NumberSharedASV")
nlsfit(data,model=4)
dta.nls <- nls(NumberSharedASV ~ (a + b * AgeMonth + c * I(AgeMonth^2)) * (AgeMonth <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (AgeMonth > -0.5 * b/c), observed,
               start=list(a=73.02278413,b=as.numeric(62.54999089),c=as.numeric(-2.13586787)))
observed$line <- predict(dta.nls)
a<-ggplot(observed,aes(x=AgeMonth, y=NumberSharedASV,colour=Weaned))+ 
  annotate("rect", xmin = 11.8, xmax = 18.2, ymin =0, ymax = 1000,alpha = 0.15,color = NA)+
  annotate("segment", x = 14.64275757, 
           xend = 14.64275757, y = 1000, 
           yend = 0,
           colour = "gray50",size=1,linetype = "dotted")+
  geom_point(size=2.5)+ 
  scale_color_manual(values=c("#440154FF", "#5DC863FF"),
                     name="",
                     breaks=c("No", "Yes"),
                     labels=c("Unweaned", "Weaned") )+
  xlab("Age (months)") + ylab("Number of shared ASVs with mothers")+
  theme(legend.position="top",legend.text =element_text(size=18,face="bold"),axis.text=element_text(size=14),axis.title=element_text(size=18,face="bold"))+
  geom_line(aes(y =line),size=1.8,color='coral')
a


data=observed[c("AgeMonth","bcval")]
nlsplot (data,model=4,xlab="Age",ylab="bcval")
nlsfit(data,model=4)
dta.nls <- nls(bcval ~ (a + b * AgeMonth + c * I(AgeMonth^2)) * (AgeMonth <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (AgeMonth > -0.5 * b/c), observed,
               start=list(a= 0.85338667,b=as.numeric( -0.05065322),c=as.numeric(0.00147669)))
observed$line <- predict(dta.nls)
b<-ggplot(observed,aes(x=AgeMonth, y=bcval,colour=Weaned))+ 
  annotate("rect", xmin = 14.4, xmax = 23.5, ymin =0.2, ymax =1,alpha = 0.15,color = NA)+
  annotate("segment", x = 17.15089360, 
           xend = 17.15089360, y = 1, 
           yend = 0.2,
           colour = "gray50",size=1,linetype = "dotted")+
  geom_point(size=2)+ 
  scale_color_manual(values=c("#440154FF", "#5DC863FF"))  +
  xlab("") + ylab("Bray-Curtis dissimilarity")+ 
  theme(legend.position="none",axis.text=element_text(size=12),axis.title=element_text(size=18,face="bold"))+
  geom_line(aes(y =line),size=1.8,color='coral')
b

data=observed[c("AgeMonth","uval")]
nlsplot (data,model=4,xlab="Age",ylab="Observed")
nlsfit(data,model=4)
dta.nls <- nls(uval ~ (a + b * AgeMonth + c * I(AgeMonth^2)) * (AgeMonth <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (AgeMonth > -0.5 * b/c), observed,
               start=list(a=0.83146173,b=as.numeric(-0.05216982),c=as.numeric(0.00170618)))
observed$line <- predict(dta.nls)
c<-ggplot(observed,aes(x=AgeMonth, y=uval,colour=Weaned))+ 
  annotate("rect", xmin = 12.6, xmax = 19.2, ymin =0.2, ymax =1,alpha = 0.15,color = NA)+
  annotate("segment", x = 15.28851965, 
           xend = 15.28851965, y = 1, 
           yend = 0.2,
           colour = "gray50",size=1,linetype = "dotted")+
  geom_point(size=2)+ 
  scale_color_manual(values=c("#440154FF", "#5DC863FF"))  +
  xlab("Age (months)") + ylab("Unweigthed Unifrac dissimilarity")+ 
  theme(legend.position="none",axis.text=element_text(size=13),axis.title=element_text(size=18,face="bold"))+
  geom_line(aes(y =line),size=1.8,color='coral')
c

data=observed[c("AgeMonth","wval")]
nlsplot (data,model=4,xlab="Age",ylab="Observed")
nlsfit(data,model=4)
dta.nls <- nls(wval ~ (a + b * AgeMonth + c * I(AgeMonth^2)) * (AgeMonth <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (AgeMonth > -0.5 * b/c), observed,
               start=list(a=1.050061e-01,b=as.numeric(-1.012423e-02),c=as.numeric(4.261800e-04)))
observed$line <- predict(dta.nls)
d<-ggplot(observed,aes(x=AgeMonth, y=wval,colour=Weaned))+
  annotate("rect", xmin = 8.9, xmax = 16.4, ymin =0.02, ymax =0.16,alpha = 0.15,color = NA)+
  annotate("segment", x = 1.187779e+01, 
           xend = 1.187779e+01, y = 0.16, 
           yend = 0.02,
           colour = "gray50",size=1,linetype = "dotted")+
  geom_point(size=2)+
  scale_color_manual(values=c("#440154FF", "#5DC863FF"))  +
  xlab("") + ylab("Weigthed Unifrac dissimilarity")+ 
  theme(legend.position="none",axis.text=element_text(size=13),axis.title=element_text(size=18,face="bold"))+
  geom_line(aes(y =line),size=1.8,color='coral')
d






######################## WHAT ARE THE SHARED ASVs in EARLY-LIFE ############################### 
load("datasets/infant and female gelada microbiomes_SILVA_ASV filter 500 reads.RData")
gelada_physeq 
otu_table(gelada_physeq)[1:5,1:5]
head(sample_data(gelada_physeq))
taxa<-data.frame(tax_table(gelada_physeq))
taxa$ASV<-rownames(taxa)
head(taxa)

#Add loading scores
pc1<-read.csv("datasets/LoadingScoresTaxa_infants.txt",sep=",",header=TRUE)
head(pc1)
taxa1<-merge(taxa,pc1[c("taxa_id","PC1")],by.x="ASV",by.y="taxa_id",all.x=TRUE,all.y=FALSE)
head(taxa1)

#Abundance and prevalence of ASV in infant samples
INF = prune_samples(sample_data(gelada_physeq)$SampleType %in% "Infant" & sample_data(gelada_physeq)$AgeMonth<=12, gelada_physeq)
INF #184
prev.i = apply(X = otu_table(INF),MARGIN = ifelse(taxa_are_rows(INF), yes = 1, no = 2),FUN = function(x){sum(x > 0)})
prev.inf = data.frame(Prevalence = prev.i, #nb de samples ou OTU est pst
                      TotalAbundance = taxa_sums(INF), #sum nb de reads per OTU belonging to the same phylum
                      tax_table(INF))
prev.inf$RA.Inf<-(prev.inf$TotalAbundance*100)/sum(prev.inf$TotalAbundance)
prev.inf$PREV.Inf<-(prev.inf$Prevalence*100)/nsamples(INF)
prev.inf$ASV<-rownames(prev.inf)
head(prev.inf)
taxa2<-merge(taxa1,prev.inf[c("ASV","RA.Inf","PREV.Inf")],all.x=TRUE,all.y=FALSE)

#Abundane and prevalence of ASV in female samples
AF= prune_samples(sample_data(gelada_physeq)$SampleType %in% "Adult female", gelada_physeq) 
AF #620
prev.a = apply(X = otu_table(AF),MARGIN = ifelse(taxa_are_rows(AF), yes = 1, no = 2),FUN = function(x){sum(x > 0)})
prev.af = data.frame(Prevalence = prev.a, #nb de samples ou OTU est pst
                     TotalAbundance = taxa_sums(AF), #sum nb de reads per OTU belonging to the same phylum
                     tax_table(AF))
prev.af$RA.AF<-(prev.af$TotalAbundance*100)/sum(prev.af$TotalAbundance)
prev.af$PREV.AF<-(prev.af$Prevalence*100)/nsamples(AF)
prev.af$ASV<-rownames(prev.af)
head(prev.af)
taxa3<-merge(taxa2,prev.af[c("ASV","RA.AF","PREV.AF")],all.x=TRUE,all.y=FALSE)
head(taxa3)


### ASV SHARED ###
df<-readRDS("datasets/mother-offspring pairs/TaxaShared_ALL PAIRS_css_ASV 500 reads.txt")
df<-droplevels(subset(df,df$AgeMonth1<=12))
length(unique(df$SampleID1))#136
df1<-data.frame(table(df$ASVShared))
colnames(df1)<-c("ASV","NbDyadsShared")
df1$PropDyadShared<-(df1$NbDyadsShared*100)/length(unique(df$SampleID1))
head(df1)

taxa4<-merge(taxa3,df1,all.x=TRUE,all.y=FALSE)
taxa4$NbDyadsShared[taxa4$NbDyadsShared %in% NA]<-0
taxa4$PropDyadShared[taxa4$PropDyadShared %in% NA]<-0
taxa5<-droplevels(subset(taxa4,taxa4$PC1 %!in% NA))
barplot(sort(taxa5$PC1,decreasing=TRUE),cex.lab=1.5,ylab="Loading scores PC1",xlab="ASVs")
abline(h=-0.025,col="red")
abline(h=0.04,col="red")

taxa5$ASV.Type<-rep("middle",length(taxa5$ASV))
taxa5$ASV.Type[taxa5$PC1 <=as.numeric(-0.025)]<-"early"
taxa5$ASV.Type[taxa5$PC1 >=as.numeric(0.04)]<-"late"
table(taxa5$ASV.Type)

taxa5$ASV.Type<-as.factor(taxa5$ASV.Type)
taxa5$Size<-rep(1,length(taxa5$ASV))
taxa5$Size[taxa5$ASV.Type %in% c("early","late")]<-3
taxa5$Shape<-rep(1,length(taxa5$ASV))
taxa5$Shape[taxa5$ASV.Type %in% c("early","late")]<-16


a<-ggplot(taxa5, aes(y=RA.Inf+0.001, x=PropDyadShared,color=ASV.Type)) + 
  xlab("Proportion of mother-infant pairs sharing the ASV (%)")+
  ylab(expression(atop(bold("Relative abundance of ASVs"), paste(bold("among infant samples (%)")))))+
  theme(legend.position="top",axis.text=element_text(size=14),
        legend.title = element_text(size=18, face="bold"),
        legend.text = element_text(size = 18),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16, face="bold"))+ 
  geom_point(alpha=0.5,size=taxa5$Size,shape=taxa5$Shape)+scale_y_log10()+
  scale_colour_manual(values=c("#440154FF","gray60","#5DC863FF"),
                      name="ASV type",
                      breaks=c("early","middle","late"),
                      labels=c("early life ASV","other","later-life ASV"))
a

b<-ggplot(taxa5, aes(y=RA.AF+0.001, x=PropDyadShared,color=ASV.Type)) + 
  xlab("Proportion of mother-infant pairs sharing the ASV (%)")+
  ylab(expression(atop(bold("Relative abundance of ASVs"), paste(bold("among adult female samples (%)")))))+
  theme(legend.position="top",axis.text=element_text(size=14),
        legend.title = element_text(size=18, face="bold"),
        legend.text = element_text(size = 18),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16, face="bold"))+ 
  geom_point(alpha=0.5,size=taxa5$Size,shape=taxa5$Shape)+scale_y_log10()+
  scale_colour_manual(values=c("#440154FF","gray70","#5DC863FF"),
                      name="ASV type",
                      breaks=c("early","middle","late"),
                      labels=c("early life ASV","other","later-life ASV"))
b


c<-ggplot(taxa5, aes(y=PREV.Inf, x=PropDyadShared,color=ASV.Type)) + 
  xlab("Proportion of mother-infant pairs sharing the ASV (%)")+
  ylab(expression(atop(bold("% infant samples"), paste(bold("containing a given ASV")))))+
  theme(legend.position="top",axis.text=element_text(size=14),
        legend.title = element_text(size=18, face="bold"),
        legend.text = element_text(size = 18),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16, face="bold"))+ 
  geom_point(alpha=0.5,size=taxa5$Size,shape=taxa5$Shape)+scale_y_log10()+
  scale_colour_manual(values=c("#440154FF","gray60","#5DC863FF"),
                      name="ASV type",
                      breaks=c("early","middle","late"),
                      labels=c("early life ASV","other","later-life ASV"))

d<-ggplot(taxa5, aes(y=PREV.AF, x=PropDyadShared,color=ASV.Type)) + 
  xlab("Proportion of mother-infant pairs sharing the ASV (%)")+
  ylab(expression(atop(bold("% adult female samples"), paste(bold("containing a given ASV")))))+
  theme(legend.position="top",axis.text=element_text(size=14),
        legend.title = element_text(size=18, face="bold"),
        legend.text = element_text(size = 18),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16, face="bold"))+ 
  geom_point(alpha=0.5,size=taxa5$Size,shape=taxa5$Shape)+scale_y_log10()+
  scale_colour_manual(values=c("#440154FF","gray60","#5DC863FF"),
                      name="ASV type",
                      breaks=c("early","middle","late"),
                      labels=c("early life ASV","other","later-life ASV"))


ggarrange(a,b,c,d,ncol=2,nrow=2,common.legend=TRUE,align="hv")





