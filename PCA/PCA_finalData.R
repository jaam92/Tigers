#Load necessary packages
library(tidyverse)
library(GENESIS)
library(GWASTools)
library(SNPRelate)

#Set working directory and load files
setwd("~/TigerProject/PCA")
#snpgdsBED2GDS(bed.fn = "highcov-only.ba-AN-MM-pcc-GM.bed", bim.fn = "highcov-only.ba-AN-MM-pcc-GM.bim", fam.fn = "highcov-only.ba-AN-MM-pcc-GM.fam", out.gdsfn = "highcov-only.ba-AN-MM-pcc-GM.gds", cvt.chr = "char")

#Open files and prep
popsDF = popsDF = read_csv("~/TigerProject/IndivFiles/individual_ids.csv") %>%
  mutate(combo = ifelse(Subspecies2 == "Unknown", paste(Subspecies2, "-", Phenotype, sep = ""), Subspecies2)) #metadata

gds = snpgdsOpen("highcov-only.ba-AN-MM-pcc-GM.gds") #open genofile
sampIds = read.gdsn(index.gdsn(gds, "sample.id")) #grab sample ids 
famIds = as.data.frame(sampIds) %>%
  left_join(popsDF, by = c("sampIds"="Individual")) %>%
  pull(Subspecies2)#make family ids

#LD prune because data set is small, set r2 = 0.7, maf 5%, and 50 snp window 
snpset = snpgdsLDpruning(gds, sample.id=sampIds, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F) 
pruned = unlist(snpset, use.names=FALSE)

#Run KING
KING = snpgdsIBDKING(gds, sample.id = sampIds, snp.id=pruned, autosome.only = F)

#Make GRM
KINGmat = KING$kinship %>%
  as.matrix()  
colnames(KINGmat) = KING$sample.id 
rownames(KINGmat) = KING$sample.id 
snpgdsClose(gds)

#Reopen GDS data 
set.seed(10192021)
genoFile = GdsGenotypeReader(filename = "highcov-only.ba-AN-MM-pcc-GM.gds")#read in GDS data
genoData = GenotypeData(genoFile)#create a GenotypeData class object

##Partition data into relateds and unrelated (less than first cousins)
sampset = pcairPartition(kinobj = KINGmat, kin.thresh=2^(-9/2), divobj = KINGmat, div.thresh=2^(-9/2)) 
unrelsID = sampset$unrels

#####Look at relatedness in full data and unrelated
# fullDataKinship = reshape2::melt(KINGmat) %>%
#   filter(Var1%in%unrelsID & Var2%in%unrelsID) %>%
#   mutate(pop1 = popsDF$Subspecies2[match(Var1, popsDF$Individual)],
#          pop2 = popsDF$Subspecies2[match(Var2, popsDF$Individual)]) %>%
#   filter(pop1 == pop2 & Var1 != Var2)  
# 
# unrelateds = reshape2::melt(KINGmat) %>%
#   filter(Var1%in%unrelsID & Var2%in%unrelsID) %>%
#   mutate(pop1 = popsDF$Subspecies2[match(Var1, popsDF$Individual)],
#          pop2 = popsDF$Subspecies2[match(Var2, popsDF$Individual)]) %>%
#   filter(pop1 == pop2 & Var1 != Var2)  
# 

####Run PC-AiR
TigerspcAir = pcair(genoData, kinobj = KINGmat, kin.thresh=2^(-9/2), divobj = KINGmat, div.thresh=2^(-9/2), unrel.set = unrelsID, snp.include = pruned, autosome.only = F)
summary(TigerspcAir)

#Prep GDS for PC-Relate
genoData = GenotypeBlockIterator(genoData, snpBlock = 20000)#take input from PC-AiR and convert to Genotype block iterator

#Run PC-Relate
TigerpcRelate = pcrelate(genoData, pcs = TigerspcAir$vectors[,1:2], training.set = TigerspcAir$unrels, ibd.probs = FALSE) #use first 2 pcs to correct kinship for population structure (aka ancestry)
pcRelateMat = pcrelateToMatrix(TigerpcRelate, scaleKin = 1) #convert pcrelate output to GRM and don't scale kinship 

#Redo PCA with ancestry adjusted kinship 
correctedTigerpcAir = pcair(genoData, kinobj= pcRelateMat, kin.thresh=2^(-9/2), 
                            divobj= KINGmat, div.thresh=2^(-9/2), 
                            snp.include = pruned, autosome.only = F) #use ancestry adjusted pcrelate GRM in the next pca
snpgdsClose(gds)
summary(correctedTigerpcAir)

##convert pcs to data frame
pcs = correctedTigerpcAir$vectors
pc.df = as.data.frame(pcs)
names(pc.df) = paste0("PC", 1:ncol(pcs))
pc.df$sample.id = row.names(pcs)
pc.df$Subspecies = popsDF$Subspecies2[match(pc.df$sample.id, popsDF$Individual)]
pc.df$pheno = popsDF$combo[match(pc.df$sample.id, popsDF$Individual)]


####Plot ancestry adjusted PCA
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Sumatran" = "cornflowerblue", "Indochinese" = "gold4", "South China" = "plum", "Unknown"="gray25", "Unknown-Orange" = "orange", "Unknown-SnowWhite" = "red", "Unknown-Golden"="gray80", "Unknown-White"="gray20")#palette

PC1vPC2 = ggplot(pc.df, aes(x=PC1, y=PC2, color=Subspecies)) + 
  geom_point(size = 4) + 
  scale_color_manual(values=cbPalette, name="Subspecies") + 
  theme_bw() + 
  labs(y=bquote('PC2' ~'('~.(round(correctedTigerpcAir$values[2], digits = 3))~'%'~')'), x=bquote('PC1'~'('~.(round(correctedTigerpcAir$values[1], digits = 3))~'%'~')')) + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))

PC2vPC3 = ggplot(pc.df, aes(x=PC2, y=PC3, color=Subspecies)) + 
  geom_point(size = 4) + 
  scale_color_manual(values=cbPalette, name="Subspecies") + 
  theme_bw() + 
  labs(y=bquote('PC3' ~'('~.(round(correctedTigerpcAir$values[3], digits = 3))~'%'~')'), x=bquote('PC2'~'('~.(round(correctedTigerpcAir$values[2], digits = 3))~'%'~')')) + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))


P1vPC2_pheno = ggplot(pc.df, aes(x=PC1, y=PC2, color=pheno)) + 
  geom_point(size = 4) + 
  scale_color_manual(values=cbPalette, name="Phenotype") + 
  theme_bw() + 
  labs(y=bquote('PC2' ~'('~.(round(correctedTigerpcAir$values[2], digits = 3))~'%'~')'), x=bquote('PC1'~'('~.(round(correctedTigerpcAir$values[1], digits = 3))~'%'~')')) + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))