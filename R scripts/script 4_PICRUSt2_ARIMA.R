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

###################################################### CHOOSE ONE DATASET ######################################################
# KEGG LEVEL 2 #
df<-read.csv("datasets/picrust2/ko_level2_relative abundance.txt",header=TRUE, row.names=1, sep=",", na.strings="NA")
df[1:5,1:5]
dim(df)

# KEGG LEVEL 3 #
#df<-read.csv("datasets/picrust2/ko_level3_relative abundance.txt",header=TRUE, row.names=1,sep=",", na.strings="NA") 
#df[1:5,1:5]
#dim(df)

# ENZYME NUMBER #
#df<-read.csv("datasets/picrust2/EC_relative abundance.txt",header=TRUE, row.names=1, sep=",", na.strings="NA")
#df<-df[,-1] #remove description
#df[1:5,1:5]
#dim(df)

#Filter out low abundant pathways (<0.01%)
length(which(rowMeans(df)>=0.01))
df_filter<-droplevels(df[which(rowMeans(df)>=0.01),])
dim(df_filter)

# z-scored the counts per pathway
stdize=function(x) {(x - mean(x))/sd(x)}
adj.data<-t(apply(t(df_filter),2,stdize)) 
adj.data[1:5,1:5]#rows=pathways, columns=samples
apply(adj.data,1,mean)#should be ~0 for all pathways
apply(adj.data,1,sd)#should be 1 for all pathways


########################################## Run ARIMA function ########################################## 
metadata<-readRDS("datasets/metadata.txt")
metadata<-droplevels(subset(metadata,metadata$SampleType %in% "Infant"))
head(metadata)

source("R scripts/arima.selection.R") #from: Márquez EJ, Chung C-H, Marches R, Rossi RJ, Nehar-Belaid D, Eroglu A, et al. Sexual-dimorphism in human immune system aging. Nat Commun. 2020;11: 751.
results <- arima.selection(adj.data = adj.data, metadata = metadata)
results$fitted.nonzero # arima output
results$predicted.nonzero # predicted values from the loess model

# Vizualization of age-related trajectories as heatmaps
arima_preds_adults <- results$predicted.nonzero 
colnames(arima_preds_adults) <- round(sort(unique(metadata$age)),1)
map<-as.matrix(arima_preds_adults,check.names = FALSE)
map[1:5,1:5]#rows=pathways
heatmap3(map, Colv = NA,cexRow=0.5,cexCol=1,showRowDendro = T, balanceColor = T,col=viridis(24), margins=c(4,12),  xlab = "Age (month)", ylab = "")  


########################################## HYERARCHICAL CLUSTERING ########################################## 
#The Elbows Method, compute total within-cluster sum of squares
fviz_nbclust(arima_preds_adults, kmeans, method = "wss", k.max = 20) + ggtitle("the Elbow Method")

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

