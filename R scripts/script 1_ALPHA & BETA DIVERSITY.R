rm(list=ls())
'%!in%'<-function(x,y)!('%in%'(x,y))
library(phyloseq)
library(tidyverse)
library(vegan)
library(plyr) 
library(reshape2)
library(ggplot2)
library(ggpubr)
library(gridExtra)
theme_set (theme_classic(base_size=15))


setwd("~/github/")

############################# DATASET #######################################
load("datasets/infant and female gelada microbiomes_SILVA_ASV filter 500 reads.RData")
gelada_physeq #3877 ASV - 1145 samples
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


### Sampling design ###
min<-sapply(unique(meta$Individual),function(x) min(meta$AgeMonth[meta$Individual==x]))
min1<-data.frame(cbind(Individual=as.character(unique(meta$Individual)),min))
min1$min<-as.numeric(as.character(min1$min))
min2<-min1[with(min1, order(-min)), ] 
meta$Individual<-factor(meta$Individual,levels = as.factor(as.character(min2$Individual)))
rt<-ggplot(meta,aes(x=AgeMonth, y=Individual,colour=Weaned,group = Individual)) +  geom_line(colour="gray60") + geom_point(size=2) + xlim(0,40)+scale_color_manual(values=c("#440154FF", "#5DC863FF"))  +
  xlab("Age (months)") + ylab("Immature ID")+ theme(axis.text.y=element_text(size=6),axis.text.x=element_text(size=15),axis.title=element_text(size=20,face="bold"))
rt1<-rt + labs(color='Immature weaned') +  theme(legend.position="top",legend.title = element_text(size=15,face="bold"),legend.text =element_text(size=15))
rt1


################################### 1. ALPHA DIVERSITY - DATASET NOT RAREFIED (for the GAMMs) ####################################

### Calculate alpha diversity indices ###
alpha<-estimate_richness(gelada_physeq,measures = c("Observed","Shannon"))

#Faith PD
library(picante)
faith<-pd(as.matrix(as.data.frame(t(otu_table(gelada_physeq)))),phy_tree(gelada_physeq),include.root=TRUE)#row=samples,columns=ASV
alpha1<-cbind(rownames(alpha),alpha,faith)
colnames(alpha1)[1]<-"SampleID"

#Add metadata
alpha2<-merge(alpha1,meta,by="SampleID",all.x=TRUE,all.y=FALSE)
alpha2$logReads<-log10(alpha2$NumberReads)
head(alpha2)
dim(alpha2)


### GAMMs ####
library(mgcv)  
library(mgcViz)

# ALL IMMATURES # 
m1<-gam(Shannon ~ Sex + Parity + MomRank  + Rain30 + MinT30 + logReads + s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re"), data=alpha2,gamma=1.4,method="REML")
summary(m1)

m2<-gam(Observed ~ Sex + Parity + MomRank + Rain30 + MinT30 + logReads+ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re"), data=alpha2,gamma=1.4,method="REML") 
summary(m2)

m3<-gam(PD ~ Sex + Parity + MomRank + Rain30 + MinT30 + logReads + s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re"), data=alpha2,gamma=1.4,method="REML") 
summary(m3)


# ONLY YOUNG IMMATURES (<12 months) # 
alpha3<-droplevels(subset(alpha2,alpha2$AgeMonth<=12))
dim(alpha3)#184 samples 

m1<-gam(Shannon ~ Sex + Parity + MomRank + Rain30 + MinT30 + logReads + s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re"), data=alpha3,gamma=1.4,method="REML")
summary(m1)

m2<-gam(Observed ~ Sex + Parity + MomRank + Rain30 + MinT30 + logReads+ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re"), data=alpha3,gamma=1.4,method="REML") 
summary(m2)

m3<-gam(PD ~ Sex + Parity + MomRank + Rain30 + MinT30 + logReads + s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re"), data=alpha3,gamma=1.4,method="REML") 
summary(m3)


# ONLY OLD INFANTS (>18 months) # 
alpha4<-droplevels(subset(alpha2,alpha2$AgeMonth>18))
dim(alpha4)#259 samples 

m1<-gam(Shannon ~ Sex + Parity + MomRank + Rain30 + MinT30 + logReads+ s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re"), data=alpha4,gamma=1.4,method="REML")
summary(m1)

m2<-gam(Observed ~ Sex + Parity + MomRank + Rain30 + MinT30 + logReads + s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re"), data=alpha4,gamma=1.4,method="REML") 
summary(m2)

m3<-gam(PD ~ Sex + Parity + MomRank + Rain30 + MinT30 + logReads + s(AgeMonth) + s(Individual,bs="re") + s(Unit,bs="re"), data=alpha4,gamma=1.4,method="REML") 
summary(m3)


# Model check #
mod<-m1 #CHANGE HERE
par(mfrow=c(2,2))  
gam.check(mod) 

viz <- getViz(mod)
plot(sm(viz, 1)) + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) + l_ciLine(mul = 5, colour = "blue", linetype = 2) + l_points(shape = 19, size = 1, alpha = 0.5) 
plot(sm(viz, 2)) + l_fitLine(colour = "red") + l_ciLine(mul = 5, colour = "blue", linetype = 2) + l_points(shape = 19, size = 1, alpha = 0.5) 
plot(sm(viz, 3)) + l_fitLine(colour = "red") + l_ciLine(mul = 5, colour = "blue", linetype = 2) + l_points(shape = 19, size = 1, alpha = 0.5) 




################################### 2. ALPHA DIVERSITY - DATASET RAREFIED (for plotting and quadratic plateau models) ####################################
#Rarefied dataset
set.seed(123)
raref <- rarefy_even_depth(gelada_physeq,sample.size=20000,replace=FALSE,trimOTUs=TRUE)

### Calculate alpha diversity indices ###
alpha.raref<-estimate_richness(raref,measures = c("Observed","Shannon"))

#Faith PD
library(picante)
faith.raref<-pd(as.matrix(as.data.frame(t(otu_table(raref)))),phy_tree(raref),include.root=TRUE)#row=samples,columns=ASV
alpha1.raref<-cbind(rownames(alpha.raref),alpha.raref,faith.raref)
head(alpha1.raref)
colnames(alpha1.raref)[1]<-"SampleID"

#Add metadata
alpha2.raref<-merge(alpha1.raref,meta,by="SampleID",all.x=TRUE,all.y=FALSE)
head(alpha2.raref)
dim(alpha2.raref)


### FIGURES ###
#AGE + WEANING STATUS
a<-ggplot(alpha2.raref,aes(x=AgeMonth, y=Observed,colour=Weaned))+ geom_point()+  xlab("Age (months)") + ylab("Observed Richness") + scale_color_manual(values=c("#440154FF", "#5DC863FF")) + theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"),legend.title = element_text(size=15,face="bold"))+ labs(color='Immature weaned') 
b<-ggplot(alpha2.raref,aes(x=AgeMonth, y=Shannon,colour=Weaned))+ geom_point()+  xlab("Age (months)") + ylab("Shannon Index") + scale_color_manual(values=c("#440154FF", "#5DC863FF")) + theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"),legend.title = element_text(size=15,face="bold"))+ labs(color='Immature weaned') 
c<-ggplot(alpha2.raref,aes(x=AgeMonth, y=PD,colour=Weaned))+ geom_point()+  xlab("Age (months)") + ylab("Faith's PD") + scale_color_manual(values=c("#440154FF", "#5DC863FF")) + theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"),legend.title = element_text(size=15,face="bold"))+ labs(color='Immature weaned') 
ggarrange(a,b,c,ncol=3,common.legend=TRUE)

#SEX
a<-ggplot(alpha2.raref,aes(x=AgeMonth, y=Observed,colour=Sex))+ geom_point() + xlab("Age (months)") + ylab("Observed Richness")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold")) +geom_smooth(method="gam",se=F)
b<-ggplot(alpha2.raref,aes(x=AgeMonth, y=Shannon,colour=Sex))+ geom_point() + xlab("Age (months)") + ylab("Shannon Index")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold")) +geom_smooth(method="gam",se=F)
c<-ggplot(alpha2.raref,aes(x=AgeMonth, y=PD,colour=Sex))+ geom_point() + xlab("Age (months)") + ylab("Faith's PD")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold")) +geom_smooth(method="gam",se=F)
ggarrange(a,b,c,ncol=3,common.legend=TRUE)

#MATERNAL PARITY
a<-ggplot(alpha2.raref,aes(x=AgeMonth, y=Observed,colour=Parity))+ geom_point() + xlab("Age (months)") + ylab("Observed Richness")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"))+geom_smooth(method="gam",se=F)
b<-ggplot(alpha2.raref,aes(x=AgeMonth, y=Shannon,colour=Parity))+ geom_point() + xlab("Age (months)") + ylab("Shannon Index")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"))+geom_smooth(method="gam",se=F)
c<-ggplot(alpha2.raref,aes(x=AgeMonth, y=PD,colour=Parity))+ geom_point() + xlab("Age (months)") + ylab("Faith's PD")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"))+geom_smooth(method="gam",se=F)
ggarrange(a,b,c,ncol=3,common.legend=TRUE)

#MATERNAL RANK
a<-ggplot(alpha2.raref,aes(x=AgeMonth, y=Observed,colour=MomRank_cat))+ geom_point() + xlab("Age (months)") + ylab("Observed Richness")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"))+geom_smooth(method="gam",se=F)
b<-ggplot(alpha2.raref,aes(x=AgeMonth, y=Shannon,colour=MomRank_cat))+ geom_point() + xlab("Age (months)") + ylab("Shannon Index")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"))+geom_smooth(method="gam",se=F)
c<-ggplot(alpha2.raref,aes(x=AgeMonth, y=PD,colour=MomRank_cat))+ geom_point() + xlab("Age (months)") + ylab("Faith's PD")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"))+geom_smooth(method="gam",se=F)
ggarrange(a,b,c,ncol=3,common.legend=TRUE)

#RAIN30
a<-ggplot(alpha2.raref,aes(x=AgeMonth, y=Observed,colour=Rain30_cat))+ geom_point() + xlab("Age (months)") + ylab("Observed Richness")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"))+geom_smooth(method="gam",se=F)
b<-ggplot(alpha2.raref,aes(x=AgeMonth, y=Shannon,colour=Rain30_cat))+ geom_point() + xlab("Age (months)") + ylab("Shannon Index")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"))+geom_smooth(method="gam",se=F)
c<-ggplot(alpha2.raref,aes(x=AgeMonth, y=PD,colour=Rain30_cat))+ geom_point() + xlab("Age (months)") + ylab("Faith's PD")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"))+geom_smooth(method="gam",se=F)
ggarrange(a,b,c,ncol=3,common.legend=TRUE)

#MINT30
a<-ggplot(alpha2.raref,aes(x=AgeMonth, y=Observed,colour=MinT30_cat))+ geom_point() + xlab("Age (months)") + ylab("Observed Richness")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"))+geom_smooth(method="gam",se=F)
b<-ggplot(alpha2.raref,aes(x=AgeMonth, y=Shannon,colour=MinT30_cat))+ geom_point() + xlab("Age (months)") + ylab("Shannon Index")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"))+geom_smooth(method="gam",se=F)
c<-ggplot(alpha2.raref,aes(x=AgeMonth, y=PD,colour=MinT30_cat))+ geom_point() + xlab("Age (months)") + ylab("Faith's PD")+theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"))+geom_smooth(method="gam",se=F)
ggarrange(a,b,c,ncol=3,common.legend=TRUE)



### QUADRATIC PLATEAU MODELS ###
#Quadratic plateau models were implemented to find the age at which alpha diversity reaches a plateau
library(easynls)

#SHANNON INDEX
data=alpha2.raref[c("AgeMonth","Shannon")]
nlsplot (data,model=4,xlab="Age",ylab="Observed")
nlsfit(data,model=4)
dta.nls <- nls(Shannon ~ (a + b * AgeMonth + c * I(AgeMonth^2)) * (AgeMonth <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (AgeMonth > -0.5 * b/c), alpha2.raref,
               start=list(a=1.9418723,b= 0.8203449,c=as.numeric(-0.0565565)))
alpha2.raref$line <- predict(dta.nls)
a1<-ggplot(alpha2.raref,aes(x=AgeMonth, y=Shannon,colour=Weaned))+ ylim(1.3,6)+
  annotate("rect", xmin = 6.4, xmax = 8.2, ymin = 1.3, ymax = 6,alpha = 0.15,color = NA)+
  annotate("segment", x = 7.2524370, 
           xend = 7.2524370, y = 6, 
           yend = 1.3,
           colour = "gray50",size=1,linetype = "dotted")+
  geom_point(size=2)+
  scale_color_manual(values=c("#440154FF", "#5DC863FF"),
                     name="",
                     breaks=c("No", "Yes"),
                     labels=c("Unweaned", "Weaned"))+
  xlab("Age (months)") + ylab("Shannon index")+ 
  theme(legend.position="top",legend.text =element_text(size=18,face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=18,face="bold"))+
  geom_line(aes(y =line),size=1.8,color='coral')
a1



#OBSERVED RICHNESS
data=alpha2.raref[c("AgeMonth","Observed")]
nlsplot (data,model=4,xlab="Age",ylab="Observed") # model=4 -> quadratic plateau model
nlsfit(data,model=4)
dta.nls <- nls(Observed ~ (a + b * AgeMonth + c * I(AgeMonth^2)) * (AgeMonth <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (AgeMonth > -0.5 * b/c), alpha2.raref,
               start=list(a= 23.85777785,b=  126.47260830,c=as.numeric(-6.19826228)))
alpha2.raref$line <- predict(dta.nls)
b1<-ggplot(alpha2.raref,aes(x=AgeMonth, y=Observed,colour=Weaned))+ylim(50,1070)+
  annotate("rect", xmin = 8.8, xmax = 11.8, ymin =50, ymax = 1070,alpha = 0.15,color = NA)+
  annotate("segment", x = 10.20226336, 
           xend = 10.20226336, y = 1070, 
           yend = 50,
           colour = "gray50",size=1,linetype = "dotted")+
  geom_point(size=2)+ 
  scale_color_manual(values=c("#440154FF", "#5DC863FF"))  +
  xlab("Age (months)") + ylab("Observed Richness")+ 
  theme(legend.position="none",axis.text=element_text(size=13),axis.title=element_text(size=18,face="bold"),plot.margin = unit(c(1,0.5,1,1), "cm"))+
  geom_line(aes(y =line),size=1.8,color='coral')
b1


#FAITH'S PD
data=alpha2.raref[c("AgeMonth","PD")]
nlsplot (data,model=4,xlab="Age",ylab="PD")#quadratic plateau
nlsfit(data,model=4)#11 months
dta.nls <- nls(PD ~ (a + b * AgeMonth + c * I(AgeMonth^2)) * (AgeMonth <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (AgeMonth > -0.5 * b/c), alpha2.raref,
               start=list(a=  6.31411533,b= 6.46010078,c=as.numeric(-0.33573014)))
alpha2.raref$line <- predict(dta.nls)
c1<-ggplot(alpha2.raref,aes(x=AgeMonth, y=PD,colour=Weaned))+ ylim(4.5,50)+
  annotate("rect", xmin = 8.5, xmax = 10.9, ymin =4.5, ymax = 50,alpha = 0.15,color = NA)+
  annotate("segment", x = 9.62097239, 
           xend = 9.62097239, y = 50, 
           yend =4.5,
           colour = "gray50",size=1,linetype = "dotted")+
  geom_point(size=2)+  
  scale_color_manual(values=c("#440154FF", "#5DC863FF"))+
  xlab("Age (months)") + ylab("Faith's PD")+ 
  theme(legend.position="none",axis.text=element_text(size=13),axis.title=element_text(size=18,face="bold"),plot.margin = unit(c(1,0.5,1,1), "cm"))+
  geom_line(aes(y =line),size=1.8,color='coral')
c1





########################################## 3. BETA DIVERSITY - DATASET NOT RAREFIED ##############################################
#Beta-diversity was computed as the Aitchison distance, which is the Euclidean distance between samples after centered log-ratio (clr) 
#transformation of the raw counts (a pseudocount is added).
ftbl<- as.matrix(t(otu_table(gelada_physeq)))
ftbl[1:5,1:5]#rows=samples
taxa<-data.frame(tax_table(gelada_physeq))
taxa$taxa_id<-rownames(taxa)

### clr-transformation of the counts ###
#https://github.com/ggloor/CoDa_microbiome_tutorial/wiki/Part-1:-Exploratory-Compositional-PCA-biplot
library(zCompositions)
czm <- cmultRepl(ftbl,label=0, method="CZM") #add a pseudo-count to the zeros (samples must be ROWS in ftbl)
czm[1:5,1:5]
clr <- t(apply(czm, 1, function(x){log(x) - mean(log(x))}))
clr[1:5,1:5] #rows=samples

### Perform the Principal Components Analysis (PCA) ###
#PCA on the Aitchison dissimilarity matrix was used to examine how immatures samples clustered by age 
pca <- prcomp(clr)
plot(pca$x[,1],pca$x[,2])

#Calculate the variance explained by PC1 and PC2
d.mvar <- sum(pca$sdev^2) # total variance
PC1 <- paste("PC1 ","(",round(sum(pca$sdev[1]^2)/d.mvar, 3)*100,"%",")",sep="")
PC2 <- paste("PC2 ","(",round(sum(pca$sdev[2]^2)/d.mvar, 3)*100,"%",")",sep="")

#Add metadata
df<-data.frame(pca$x)
df$SampleID<-rownames(df)
df1<-merge(df,meta,by="SampleID",all.x=TRUE,all.y=FALSE)
head(df1)


### Loading scores of ASVs ###
pca$rotation[1:5,1:5]
t1<-data.frame(pca$rotation)
t1$taxa_id<-rownames(t1)
barplot(sort(t1$PC1,decreasing=TRUE),cex.lab=1.5,ylab="Loading scores PC1",xlab="ASVs")
taxa1<-merge(taxa,t1[c("taxa_id","PC1")],by="taxa_id")
taxa2<-taxa1[with(taxa1,order(-PC1)), ] 
head(taxa2)


### PERMANOVA TESTS (takes a while) ###
aich<-dist(clr) #Aitchison distance
meta$log10_NumberReads<-log10(meta$NumberReads)
#all_null<-adonis2(aich ~ Individual + log10_NumberReads, data = meta, permutations = 10000, by = "margin")
perm <- how(nperm = 10000)
setBlocks(perm) <- with(meta,Individual)
#all1<-adonis2(aich ~ log10_NumberReads + Unit + AgeMonth + Sex + Parity + MomRank + Rain30 + MinT30, data = meta, permutations = perm, by = "margin")
#all2<-adonis2(aich ~ log10_NumberReads + Unit + AgeMonth * Parity + AgeMonth * MomRank + Sex + Rain30 + MinT30, data = meta, permutations = perm, by = "margin")



### PLOT BETA DIVERSITY & AGE ###
#AGE CATEGORICAL
a<-ggplot(df1,aes(x=PC1,y=PC2,colour=AgeMonth_cat)) + geom_point(size=1.5)+xlab(PC1) + ylab(PC2)+
  scale_color_manual(values=c("#440154FF","#D95F02","#E7298A","#E69F00", "#1F78B4","black", "#5DC863FF"))+
  theme(legend.position="right",axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=12))
a1<-a + guides(color=guide_legend("Age (months)")) +  theme(legend.title = element_text(size=10,face="bold"),legend.text =element_text(size=10))
a1

#AGE CONTINUOUS
a<-ggplot(df1,aes(y=PC1, x=AgeMonth,colour=Weaned))+ geom_point()+  xlab("Age (months)") + ylab(PC1) + scale_color_manual(values=c("#440154FF", "#5DC863FF")) + theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"),legend.title = element_text(size=15,face="bold"))+ labs(color='Immature weaned') 
b<-ggplot(df1,aes(y=PC2, x=AgeMonth,colour=Weaned))+ geom_point()+  xlab("Age (months)") + ylab(PC2) + scale_color_manual(values=c("#440154FF", "#5DC863FF")) + theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"),legend.title = element_text(size=15,face="bold"))+ labs(color='Immature weaned') 
ggarrange(a,b,ncol=2,common.legend=TRUE)
cor.test(df1$PC1,df1$AgeMonth, alternative = "two.sided", method = "pearson")


### QUADRATIC PLATEAU MODELS - PC1 & AGE ###
#A quadratic plateau model was implemented to find the age at which Aitchison beta diversity reaches a plateau
data=df1[c("AgeMonth","PC1")]
nlsplot (data,model=4,xlab="Age",ylab="PC1")
nlsfit(data,model=4)
dta.nls <- nls(PC1 ~ (a + b * AgeMonth + c * I(AgeMonth^2)) * (AgeMonth <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (AgeMonth > -0.5 * b/c), df1,
               start=list(a=as.numeric(-102.8260794),b=  14.4411188,c=as.numeric( -0.4190814)))
df1$line <- predict(dta.nls)
d1<-ggplot(df1,aes(x=AgeMonth, y=PC1,colour=Weaned))+ 
  annotate("rect", xmin = 15.5, xmax = 19.4, ymin =-116, ymax = 100,alpha = 0.15,color = NA)+
  annotate("segment", x = 17.2294921, 
           xend = 17.2294921, y = 100, 
           yend = -116,
           colour = "gray50",size=1,linetype = "dotted")+
  geom_point(size=2)+ 
  scale_color_manual(values=c("#440154FF", "#5DC863FF"),
                     name="",
                     breaks=c("No", "Yes"),
                     labels=c("Unweaned", "Weaned") )+
  xlab("Age (months)") + ylab(PC1)+ 
  theme(legend.position="top",legend.text =element_text(size=18,face="bold"),axis.text=element_text(size=12),axis.title=element_text(size=18,face="bold"))+
  geom_line(aes(y =line),size=1.8,color='coral')
d1



### BETA DIVERSITY & OTHER COVARIATES ###
#SEQUENCING DEPTH
ggplot(df1,aes(y=PC1, x=NumberReads))+ geom_point() + ylab(PC1) + xlab("Sequencing depth") + scale_x_log10()
ggplot(df1,aes(y=PC2, x=NumberReads))+ geom_point() + ylab(PC2) + xlab("Sequencing depth") + scale_x_log10()

#SEX
ggplot(df1,aes(x=PC1,y=PC2,colour=Sex)) + geom_point() + xlab(PC1) + ylab(PC2)
ggplot(df1, aes(y=PC1,x=AgeMonth,colour=Sex)) + geom_point() + ylab(PC1) + geom_smooth(method="gam",se=F)

#MATERNAL PARITY
ggplot(df1,aes(x=PC1,y=PC2,colour=Parity)) + geom_point() + xlab(PC1) + ylab(PC2)
ggplot(df1, aes(y=PC1,x=AgeMonth,colour=Parity)) + geom_point() + ylab(PC1) +geom_smooth(method="gam",se=F)

#MATERNAL RANK
ggplot(df1,aes(x=PC1,y=PC2,colour=MomRank_cat)) + geom_point() + xlab(PC1) + ylab(PC2)
ggplot(df1, aes(y=PC1,x=AgeMonth,colour=MomRank_cat)) + geom_point() + ylab(PC1) +geom_smooth(method="gam",se=F)

#RAIN30
ggplot(df1,aes(x=PC1,y=PC2,colour=Rain30_cat)) + geom_point() + xlab(PC1) + ylab(PC2)
ggplot(df1, aes(y=PC1,x=AgeMonth,colour=Rain30_cat)) + geom_point() + ylab(PC1) +geom_smooth(method="gam",se=F)

#MINT30
ggplot(df1,aes(x=PC1,y=PC2,colour=MinT30_cat)) + geom_point() + xlab(PC1) + ylab(PC2)
ggplot(df1, aes(y=PC1,x=AgeMonth,colour=MinT30_cat)) + geom_point() + ylab(PC1) +geom_smooth(method="gam",se=F)




########################################## 4. BETA DIVERSITY - CSS  TRANSFORMATION ##############################################
#We replicated those PERMANOVA analyses using more classical measures of beta diversity (unweighted and weighted UniFrac dissimilarity) on a dataset where counts have been CSS transformed

### CSS TRANSFORNATION ###
library(metagenomeSeq)
otu_table(gelada_physeq)[1:5,1:5]
OTU_read_count = as.data.frame(otu_table(gelada_physeq))#rows=taxa and columns=samples

#Convert OTU table into package format
metaSeqObject= newMRexperiment(OTU_read_count) 

#CSS normalization
metaSeqObject_CSS  = cumNorm(metaSeqObject,p=cumNormStatFast(metaSeqObject))

#convert CSS normalized data into data.frame-formatted OTU table (log transformed data)
OTU_read_count_CSS = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE))
OTU_read_count_CSS[1:5,1:5]

css<-gelada_physeq
otu_table(css)<-otu_table(OTU_read_count_CSS,taxa_are_rows=TRUE)
otu_table(css)[1:5,1:5]

meta<-data.frame(sample_data(css))
meta$log10_NumberReads<-log10(meta$NumberReads)


### Bray-Curtis dissimilarity ###
bcdist = phyloseq::distance(css, method="bray",normalized=TRUE, parallel=FALSE, fast=TRUE)  
ord1 <- ordinate(css, method = "PCoA", distance = bcdist)
meta$bc.PC1=ord1$vectors[,1]
meta$bc.PC2=ord1$vectors[,2]
plot_ordination(css, ord1, color="AgeMonth_cat",axes=1:2) + geom_point(size = 1) 
a<-ggplot(meta,aes(y=bc.PC1, x=AgeMonth,colour=Weaned))+ geom_point()+  xlab("Age (months)")  + scale_color_manual(values=c("#440154FF", "#5DC863FF")) + theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"),legend.title = element_text(size=15,face="bold"))+ labs(color='Immature weaned') 
b<-ggplot(meta,aes(y=bc.PC2, x=AgeMonth,colour=Weaned))+ geom_point()+  xlab("Age (months)")  + scale_color_manual(values=c("#440154FF", "#5DC863FF")) + theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"),legend.title = element_text(size=15,face="bold"))+ labs(color='Immature weaned') 
ggarrange(a,b,ncol=2,common.legend=TRUE)
#mod1_null<-adonis2(bcdist ~ Individual, data = meta, permutations = 10000, by="margin")
#perm <- how(nperm = 10000)
#setBlocks(perm) <- with(meta,Individual)
#mod1a<-adonis2(bcdist ~ Unit + AgeMonth + Sex + Parity + MomRank + Rain30 + MinT30, data = meta, permutations = perm, by = "margin")
#mod1b<-adonis2(bcdist ~ Unit + AgeMonth * Parity + AgeMonth * MomRank + Sex + Rain30 + MinT30, data = meta, permutations = perm, by = "margin")


### Unweighted UniFrac dissimilarity ###
udist=UniFrac(css,weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)
ord2 <- ordinate(css, method = "PCoA", distance = udist)
meta$uu.PC1=ord2$vectors[,1]
meta$uu.PC2=ord2$vectors[,2]
plot_ordination(css, ord2, color="AgeMonth_cat",axes=1:2) + geom_point(size = 1) 
a<-ggplot(meta,aes(y=uu.PC1, x=AgeMonth,colour=Weaned))+ geom_point()+  xlab("Age (months)")  + scale_color_manual(values=c("#440154FF", "#5DC863FF")) + theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"),legend.title = element_text(size=15,face="bold"))+ labs(color='Immature weaned') 
b<-ggplot(meta,aes(y=uu.PC2, x=AgeMonth,colour=Weaned))+ geom_point()+  xlab("Age (months)")  + scale_color_manual(values=c("#440154FF", "#5DC863FF")) + theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"),legend.title = element_text(size=15,face="bold"))+ labs(color='Immature weaned') 
ggarrange(a,b,ncol=2,common.legend=TRUE)
#mod2_null<-adonis2(udist ~ Individual, data = meta, permutations = 10000, by="margin")
#perm <- how(nperm = 10000)
#setBlocks(perm) <- with(meta,Individual)
#mod2a<-adonis2(udist ~ Unit + AgeMonth + Sex + Parity + MomRank + Rain30 + MinT30, data = meta, permutations = perm, by = "margin")
#mod2b<-adonis2(udist ~ Unit + AgeMonth * Parity + AgeMonth * MomRank + Sex + Rain30 + MinT30, data = meta, permutations = perm, by = "margin")


### Weighted UniFrac dissimilarity ###
wdist=UniFrac(css, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
ord3 <- ordinate(css, method = "PCoA", distance = wdist)
meta$wu.PC1=ord3$vectors[,1]
meta$wu.PC2=ord3$vectors[,2]
plot_ordination(css, ord3, color="AgeMonth_cat",axes=1:2) + geom_point(size = 1) 
a<-ggplot(meta,aes(y=wu.PC1, x=AgeMonth,colour=Weaned))+ geom_point()+  xlab("Age (months)")  + scale_color_manual(values=c("#440154FF", "#5DC863FF")) + theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"),legend.title = element_text(size=15,face="bold"))+ labs(color='Immature weaned') 
b<-ggplot(meta,aes(y=wu.PC2, x=AgeMonth,colour=Weaned))+ geom_point()+  xlab("Age (months)")  + scale_color_manual(values=c("#440154FF", "#5DC863FF")) + theme(axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold"),legend.title = element_text(size=15,face="bold"))+ labs(color='Immature weaned') 
ggarrange(a,b,ncol=2,common.legend=TRUE)
#mod3_null<-adonis2(wdist ~ Individual, data = meta, permutations = 10000, by="margin")
#perm <- how(nperm = 10000)
#setBlocks(perm) <- with(meta,Individual)
#mod3a<-adonis2(wdist ~ Unit + AgeMonth + Sex + Parity + MomRank + Rain30 + MinT30, data = meta, permutations = perm, by = "margin")
#mod3b<-adonis2(wdist ~ Unit + AgeMonth * Parity + AgeMonth * MomRank + Sex + Rain30 + MinT30, data = meta, permutations = perm, by = "margin")




########################################## 5. UNIQUENESS - CSS DATASET ##############################################
### Uniqueness Weighted UniFrac distance ###
wdist=UniFrac(css, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
m3=as.matrix(wdist)
m3[upper.tri(m3,diag=T)] <-NA
jd3<-melt(m3)
colnames(jd3)<-c("Ind1","Ind2","value")
jd3<-droplevels(subset(jd3,jd3$value %!in% NA))
Min<-sapply(unique(jd3$Ind1),function(x) min(jd3$value[jd3$Ind1==x]))
jd3.1<-data.frame(cbind(SampleID=as.character(unique(jd3$Ind1)),Min))
jd3.1$Min<-as.numeric(as.character(jd3.1$Min))
jd3.2<-merge(jd3.1,meta[c("SampleID","Individual","AgeMonth","Weaned")],by="SampleID",all.x=TRUE,all.y=FALSE)

#Quadratic plateau model
data=jd3.2[c("AgeMonth","Min")]
nlsplot (data,model=4,xlab="Age",ylab="NumberSharedASV")
nlsfit(data,model=4)
dta.nls <- nls(Min ~ (a + b * AgeMonth + c * I(AgeMonth^2)) * (AgeMonth <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (AgeMonth > -0.5 * b/c),jd3.2,
               start=list(a=5.405087e-02,b=as.numeric(-2.397320e-03),c=as.numeric(6.642000e-05)))
jd3.2$line <- predict(dta.nls)
u3<-ggplot(jd3.2,aes(x=AgeMonth, y=Min,colour=Weaned))+ ylim(0.01,0.085)+
  annotate("rect", xmin = 14.9, xmax = 22.5, ymin =0.01, ymax = 0.085,alpha = 0.15,color = NA)+
  annotate("segment", x = 1.804658e+01, 
           xend = 1.804658e+01, y = 0.085, 
           yend = 0.01,
           colour = "gray50",size=1,linetype = "dotted")+
  geom_point(size=2.5)+
  scale_color_manual(values=c("#440154FF", "#5DC863FF"),
                     name="",
                     breaks=c("No", "Yes"),
                     labels=c("Unweaned", "Weaned")) +
  xlab("Age (months)") + ylab("Weighted UniFrac uniqueness")+ 
  theme(legend.position="top",legend.text =element_text(size=18,face="bold"),axis.text=element_text(size=14),axis.title=element_text(size=18,face="bold"))+
  geom_line(aes(y =line),size=1.8,color='coral')
u3



