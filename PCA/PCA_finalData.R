#Load necessary packages
library(tidyverse)
library(GENESIS)
library(GWASTools)
library(SNPRelate)

#Set working directory and load files
setwd("~/Documents/Tigers/PCA/plinkFiles/")
set.seed(10192021)
snpgdsBED2GDS(bed.fn = "highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.bed", bim.fn = "highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.bim", fam.fn = "highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.fam", out.gdsfn = "highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.gds", cvt.chr = "char")

#Open files and prep
popsDF = read_delim("~/Documents/Tigers/IndivFiles/TableX-SampleDetails.txt") 

gds = snpgdsOpen("highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.gds") #open genofile
sampIds = read.gdsn(index.gdsn(gds, "sample.id")) #grab sample ids 
famIds = as.data.frame(sampIds) %>%
  left_join(popsDF, by = c("sampIds"="Sample")) %>%
  pull(Subspecies_GroupID_Corrected)#make family ids

#LD prune because data set is small, set r2 = 0.7, maf 5%, and 50 snp window 
snpset = snpgdsLDpruning(gds, sample.id=sampIds, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F) 
pruned = unlist(snpset, use.names=FALSE)

#Run KING
KING = snpgdsIBDMLE(gds, sample.id = sampIds, snp.id=pruned, autosome.only = F)

#Make GRM
KINGmat = KING$kinship %>%
  as.matrix()  
colnames(KINGmat) = KING$sample.id 
rownames(KINGmat) = KING$sample.id 
snpgdsClose(gds)

#Reopen GDS data 
genoFile = GdsGenotypeReader(filename = "highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.updatedFID.gds")#read in GDS data
genoData = GenotypeData(genoFile)#create a GenotypeData class object

##Partition data into relateds and unrelated (less than first cousins)
sampset = pcairPartition(kinobj = KINGmat, kin.thresh=2^(-9/2), divobj = KINGmat, div.thresh=2^(-9/2)) 
unrelsID = sampset$unrels

#Output unrelateds
# x = cbind.data.frame(unrelsID) %>%
#   mutate(popsDF$Subspecies_GroupID_Corrected[match(unrelsID, popsDF$Sample)],
#          popsDF$`Depth (post filtering)`[match(unrelsID, popsDF$Sample)]) 
# names(x) = c("Sample", "Population", "Depth")
# 
# n10 = x %>% 
#   filter(Depth >= 5) %>% 
#   group_by(Population) %>% 
#   top_n(10)
# 
# n6 = n10 %>% 
#   group_by(Population) %>% 
#   top_n(6)
# 
# x$N10 = ifelse(x$Sample%in%n10$Sample, 1, 0)
# x$N6 = ifelse(x$Sample%in%n6$Sample, 1, 0)
# 
# write.table(x, file = "~/Documents/Tigers/IndivFiles/N10-N6_unrelateds.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#####Look at relatedness in full data and unrelated
# fullDataKinship = reshape2::melt(KINGmat) %>%
#   filter(Var1%in%unrelsID & Var2%in%unrelsID) %>%
#   mutate(pop1 = popsDF$Subspecies_GroupID_Corrected[match(Var1, popsDF$Sample)],
#          pop2 = popsDF$Subspecies_GroupID_Corrected[match(Var2, popsDF$Sample)]) %>%
#   filter(pop1 == pop2 & Var1 != Var2)
# 
# unrelateds = reshape2::melt(KINGmat) %>%
#   filter(Var1%in%unrelsID & Var2%in%unrelsID) %>%
#   mutate(pop1 = popsDF$Subspecies_GroupID_Corrected[match(Var1, popsDF$Sample)],
#          pop2 = popsDF$Subspecies_GroupID_Corrected[match(Var2, popsDF$Sample)]) %>%
#   filter(pop1 == pop2 & Var1 != Var2)
 

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
pcs = TigerspcAir$vectors
pc.df = as.data.frame(pcs)
names(pc.df) = paste0("PC", 1:ncol(pcs))
pc.df$sample.id = row.names(pcs)
pc.df$Subspecies = popsDF$Subspecies2[match(pc.df$sample.id, popsDF$Individual)]
pc.df$pheno = popsDF$combo[match(pc.df$sample.id, popsDF$Individual)]


####Plot ancestry adjusted PCA
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Sumatran" = "cornflowerblue", "Indochinese" = "gold4", "South China" = "plum", "Generic"="gray25", "Generic-Orange" = "orange", "Generic-SnowWhite" = "red", "Generic-Golden"="gray80", "Generic-White"="gray20")#palette

PC1vPC2 = ggplot(pc.df, aes(x=PC1, y=PC2, color=Subspecies)) + 
  geom_point(size = 4) + 
  scale_color_manual(values=cbPalette, name="Subspecies") + 
  theme_bw() + 
  labs(y=bquote('PC2' ~'('~.(round(TigerspcAir$values[2], digits = 3))~'%'~')'), x=bquote('PC1'~'('~.(round(TigerspcAir$values[1], digits = 3))~'%'~')')) + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))

PC2vPC3 = ggplot(pc.df, aes(x=PC2, y=PC3, color=Subspecies)) + 
  geom_point(size = 4) + 
  scale_color_manual(values=cbPalette, name="Subspecies") + 
  theme_bw() + 
  labs(y=bquote('PC3' ~'('~.(round(TigerspcAir$values[3], digits = 3))~'%'~')'), x=bquote('PC2'~'('~.(round(TigerspcAir$values[2], digits = 3))~'%'~')')) + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))


PC1vPC2_pheno = ggplot(pc.df, aes(x=PC1, y=PC2, color=pheno)) + 
  geom_point(size = 4) + 
  scale_color_manual(values=cbPalette, name="Phenotype") + 
  theme_bw() + 
  labs(y=bquote('PC2' ~'('~.(round(TigerspcAir$values[2], digits = 3))~'%'~')'), x=bquote('PC1'~'('~.(round(TigerspcAir$values[1], digits = 3))~'%'~')')) + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
