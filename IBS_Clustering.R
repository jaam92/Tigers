#Load necessary packages
library(tidyverse)
library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(corrplot)
library(ggdendro)
library(dendextend)
library(seriation) #olo function
library(cluster) #use silhouette

#Set working directory and load files
setwd("~/Documents/Tigers/PCA/")
set.seed(10192021)

#Open files and prep
popsDF = read_delim("~/Documents/Tigers/IndivFiles/TableX-SampleDetails.txt") 

gds = snpgdsOpen("~/Documents/Tigers/PCA/plinkFiles/highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.gds") #open genofile
sampIds = read.gdsn(index.gdsn(gds, "sample.id")) #grab sample ids 

famIds = as.data.frame(sampIds) %>%
  left_join(popsDF, by = c("sampIds"="Sample")) %>%
  pull(Subspecies_GroupID_Corrected)#make family ids

#Run just the wild tiger data
gdsSampIDs = read.gdsn(index.gdsn(gds, "sample.id"))
captive = popsDF %>%
  filter(Subspecies_GroupID_Corrected == "Generic") %>%
  pull(Sample)
samp.id.captive = unlist(gdsSampIDs[gdsSampIDs%in%captive])



#LD prune because data set is small, set r2 = 0.7, maf 5%, and 50 snp window 
snpset = snpgdsLDpruning(gds, sample.id=samp.id.captive, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F) 
pruned.captive = unlist(snpset, use.names=FALSE)


#Run KING
KING = snpgdsIBDKING(gds, sample.id = samp.id.captive, snp.id=pruned.captive, autosome.only = F, num.thread = 4)

#Make the IBS0 (probability IBS=0) matrix
kingIBSMatrix = KING$IBS0
colnames(kingIBSMatrix) <- rownames(kingIBSMatrix) <- KING$sample.id
disMatrix = as.dist(kingIBSMatrix) #convert to distance matrix

#Perform heirarchical clustering
hc_avg = hclust(disMatrix, method = "average")
hc_avg_reorder = reorder(hc_avg, disMatrix, method = "OLO") #reorder with optimal leaf ordering

#Function to plot 
dend = hc_avg_reorder %>% 
  as.dendrogram %>%
  #set("branches_k_color", k = 1, c("#0072B2", "#882255", "darkseagreen4", "gray25", "gold4", "#D55E00")) %>% 
  set("branches_lwd", 0.7) %>%
  set("labels_cex", 1.5) %>% 
  #set("labels_colors", k = 1, c("#0072B2", "#882255", "darkseagreen4", "gray25", "gold4", "#D55E00")) %>%
  set("leaves_pch", 19) %>% 
  set("leaves_cex", 0.5) 

ggd1 = as.ggdend(dend)

IBSTree = ggplot(ggd1, horiz = F) + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 30),
        axis.title=element_text(size = 30))

#SNPRelate IBS 
snpgdsIBS = snpgdsIBS(gds, sample.id = samp.id.captive, snp.id=pruned.captive, autosome.only = F, num.thread=4)
snpgdsIBSMatrix = snpgdsIBS$ibs
colnames(snpgdsIBSMatrix) <- rownames(snpgdsIBSMatrix) <- snpgdsIBS$sample.id
ibs.hc = snpgdsHCluster(snpgdsIBS)
snpgdsClose(gds)
 
 
#IBS distance matrix from snpgds 
rv = snpgdsCutTree(ibs.hc)
IBSTreeSNPgds = ggplot(rv$dendrogram, horiz = F) + 
ggtitle("IBS Distance Matrix from SNP Relate") + 
  theme(axis.text.x = element_blank(), 
         axis.text.y = element_text(size  = 20),
         axis.title=element_text(size=24))


#Check whether cluster membership is well supported (Sw ~ 0.6, so it's pretty good)
####For interpretation of silhouette width https://www.stat.berkeley.edu/~spector/s133/Clus.html
####Code: https://www.rdocumentation.org/packages/cluster/versions/2.0.8/topics/silhouette

c6 = c("tomato", "darkseagreen4", "steelblue", "purple2", "goldenrod4", "gray25")
par(mfrow= c(3,2), oma= c(0,0, 3, 0),mgp= c(1.6,.8,0), mar= .1+c(4,2,2,2))
for(k in 2:6){
plot(silhouette(cutree(hc_avg_reorder, k= k),disMatrix), 
     main = paste("k = ",k), 
     do.n.k=FALSE, 
     col = c6[1:k])
}

# Visualize silhouhette information
require("cluster")
sil <- silhouette(cutree(hc_avg_reorder, k=6),disMatrix)
plot(sil, cex.names = par("cex.axis"))


###All individuals
setwd("~/Documents/Tigers/PCA/")

gds = snpgdsOpen("high-corr-nodups-biallelic-AN-MM-pcc.gds") #open genofile
sampIds = read.gdsn(index.gdsn(gds, "sample.id")) #grab sample ids 
famIds = as.data.frame(sampIds) %>%
  left_join(popsDF, by = c("sampIds"="Individual")) %>%
  pull(Subspecies2)#make family ids

#LD prune because data set is small, set r2 = 0.7, maf 5%, and 50 snp window 
snpset = snpgdsLDpruning(gds, sample.id=sampIds, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F) 
pruned = unlist(snpset, use.names=FALSE)

#Run KING
KING = snpgdsIBDKING(gds, sample.id = sampIds, snp.id=pruned, autosome.only = F)

#Make the IBS0 (probability IBS=0) matrix
kingIBSMatrix = KING$IBS0
colnames(kingIBSMatrix) <- rownames(kingIBSMatrix) <- KING$sample.id
disMatrix = as.dist(kingIBSMatrix) #convert to distance matrix


#SNPRelate IBS 
snpgdsIBS = snpgdsIBS(gds, sample.id = sampIds, snp.id=pruned, num.thread=4)
snpgdsIBSMatrix = snpgdsIBS$ibs
colnames(snpgdsIBSMatrix) <- rownames(snpgdsIBSMatrix) <- snpgdsIBS$sample.id
ibs.hc = snpgdsHCluster(snpgdsIBS)
snpgdsClose(gds)