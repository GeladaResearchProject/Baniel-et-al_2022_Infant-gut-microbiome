rm(list=ls())
'%!in%'<-function(x,y)!('%in%'(x,y))
library(phyloseq)
library(zoo)
library(forecast)
library(reshape2)
library(parallel)
library(dplyr)
library(tidyverse)
library(viridis)
library(factoextra)
library(ggpubr)
library(heatmap3)
library(ggplot2)
library(gridExtra)
theme_set (theme_classic(base_size=15))

setwd("~/github/")

############################# DATASET #######################################
load("datasets/infant and female gelada microbiomes_SILVA_ASV filter.RData")
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


#Choose taxonomic level of analysis & aggregate the counts
level="Family"
#level="Genus"
glom = tax_glom(gelada_physeq,level,NArm=FALSE) #keep all taxa (even with NA) here for clr transformation


### clr-transformation of the counts ###
#https://github.com/ggloor/CoDa_microbiome_tutorial/wiki/Part-1:-Exploratory-Compositional-PCA-biplot
library(zCompositions)
czm <- cmultRepl(ftbl,label=0, method="CZM") #add a pseudo-count to the zeros (samples must be ROWS in ftbl)
czm[1:5,1:5]
clr <- t(apply(czm, 1, function(x){log(x) - mean(log(x))}))
clr[1:5,1:5] #rows=samples

#Filter low abundant taxa (<0.01%)
glom_ra = transform_sample_counts(glom, function(x) (x*100) / sum(x)) #transform to relative abundance
glom_ra = filter_taxa(glom_ra, function(x) mean(x)>=0.01, TRUE) 
glom_ra
taxa1<-data.frame(tax_table(glom_ra))
taxa1$Taxa<-taxa1[ ,which(colnames(taxa1) %in% level)]
taxa2<- droplevels(subset(taxa1,taxa1$Taxa %!in% NA))
dim(taxa2)
clr1<-clr[,which(colnames(clr) %in% rownames(taxa2))]
dim(clr1)
clr1[1:5,1:5]

# z-scored the counts per taxon
stdize=function(x) {(x - mean(x))/sd(x)}
adj.data<-t(apply(clr1,2,stdize)) 
adj.data[1:5,1:5]#rows=taxa
rownames(adj.data)<-taxa2$Taxa[match(rownames(adj.data),rownames(taxa2))]
dim(adj.data)
apply(adj.data,1,mean)#should be ~0 for all taxa
apply(adj.data,1,sd)#should be 1 for all taxa


########################################## Run ARIMA function ########################################## 
metadata<-readRDS("datasets/metadata.txt")
metadata<-droplevels(subset(metadata,metadata$SampleType %in% "Infant"))
head(metadata)

source("R scripts/arima.selection.R") #from: Márquez EJ, Chung C-H, Marches R, Rossi RJ, Nehar-Belaid D, Eroglu A, et al. Sexual-dimorphism in human immune system aging. Nat Commun. 2020;11: 751.
results <- arima.selection(adj.data = adj.data , metadata = metadata)
results$fitted.nonzero [1:5,1:5] # arima output
results$predicted.nonzero [1:5,1:5] # the predicted values from the loess model

# Vizualization of age-related trajectories as heatmaps
arima_preds_adults <- results$predicted.nonzero 
dim(arima_preds_adults)
colnames(arima_preds_adults) <- round(sort(unique(metadata$age)),1)
map<-as.matrix(arima_preds_adults,check.names = FALSE)
map[1:5,1:5]#rows=taxa
heatmap3(map, Colv = NA,cexRow=0.5,cexCol=1,showRowDendro = T, balanceColor = T,col=viridis(24), margins=c(4,12),  xlab = "Age (month)", ylab = "")  


########################################## HYERARCHICAL CLUSTERING ########################################## 
#The Elbows Method, compute total within-cluster sum of squares
fviz_nbclust(arima_preds_adults, kmeans, method = "wss", k.max = 20) + theme_minimal() + ggtitle("the Elbow Method")

# Cluster age-related loess trajectories
clus.distance.arima = as.dist(1 - cor(t(arima_preds_adults), use = "pa")) # calc. distances using heatmap3 formula 
h.clust.arima = hclust(clus.distance.arima,method='complete') # calc. clusters 
plot(h.clust.arima , cex = 0.6) 
rect.hclust(h.clust.arima, k = 4, border = 2:6) 

# Separate age-trajectories into k clusters 
k=4
hc.clusters.arima= data.frame(gene=rownames(arima_preds_adults),cluster=cutree(h.clust.arima,k=k),stringsAsFactors=FALSE)
table(hc.clusters.arima$cluster)

p<- list()
for (i in 1:k){
  c1<-rownames(subset(hc.clusters.arima ,cluster == i))
  mat1<-as.matrix(arima_preds_adults[rownames(arima_preds_adults) %in% c1,])
  cluster1<-melt(mat1)
  colnames(cluster1)<-c("Taxa","Age","LOESS_Pred")
  p[[i]]<-ggplot(cluster1, aes(x = Age, y = LOESS_Pred, group = Taxa, color = Taxa)) + 
    geom_line() + ylab("Predicted clr-abundance") + xlab("Age (months)") + theme(axis.title= element_text(face="bold"),legend.position = "none")
}
my_layout <- rbind(c(1:2),c(3:4))
grid.arrange(grobs=p[1:4],layout_matrix = my_layout)  
