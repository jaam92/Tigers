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
setwd("~/Tigers/PCA/plinkFiles")
set.seed(10192021)

#Open files
ancestryInfo = read_csv("~/Tigers/IndivFiles/plotting_metadata.csv")
popsDF = read_delim("~/Tigers/IndivFiles/TableX-SampleDetails.txt") 
gds = snpgdsOpen("~/Tigers/PCA/plinkFiles/highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.gds") #open genofile

#Get sample and family ids
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

#SNPRelate IBS 
snpgdsIBS = snpgdsIBS(gds, sample.id = samp.id.captive, snp.id=pruned.captive, autosome.only = F, num.thread=4)
snpgdsIBSMatrix = snpgdsIBS$ibs
colnames(snpgdsIBSMatrix) <- rownames(snpgdsIBSMatrix) <- snpgdsIBS$sample.id
ibs.hc = snpgdsHCluster(snpgdsIBS) #cluster based on IBS
snpgdsClose(gds) #close file
 
#Plot as a dendrogram
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Indochinese" = "gold4", "Sumatran" = "cornflowerblue", "South China" = "plum", "Generic"="gray25")#palette

dend = dendro_data(ibs.hc$dendrogram) 
labs = bind_cols(filter(segment(dend), x == xend & x%%1 == 0), label(dend)) %>%
  select(xend, yend, label) %>%
  mutate(ancestry = ancestryInfo$Top_Ancestry[match(label, ancestryInfo$Individual)]) 
  
#plot it 
ggplot(segment(dend)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_point(data=labs,
             aes(x=xend, y=yend, colour=ancestry), size = 3) +
  scale_colour_manual(values=cbPalette) +
  coord_flip() + 
  scale_y_reverse() +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.x=element_text(size = 20),
        axis.title.x=element_text(size = 24),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())


###All individuals
setwd("~/Tigers/PCA/")

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
