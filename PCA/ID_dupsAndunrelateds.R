#Load necessary packages
library(tidyverse)
library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(data.table)

#Set working directory and load files
setwd("~/Documents/Tigers/PCA")

#Open files and prep
popsDF = popsDF = read_csv("~/Documents/Tigers/IndivFiles/individual_ids.csv") %>%
  mutate(combo = ifelse(Subspecies2 == "Generic", paste(Subspecies2, "-", Phenotype, sep = ""), Subspecies2)) #metadata

#####Make gds files
#files <- list.files(path="~/Documents/Tigers/PCA/plinkFiles", pattern="*.bed", full.names=TRUE, recursive=FALSE)
# 
#for (f in files) {
#   print(f)
#   
#   inFile = gsub(".bed", "", f)
#   
#   snpgdsBED2GDS(bed.fn = f, bim.fn = paste0(inFile, ".bim"), fam.fn = paste0(inFile, ".fam"), out.gdsfn = paste0(inFile, ".gds"), cvt.chr = "char")
#   
# }

#####Run operations on gds files
#files <- list.files(path="~/Documents/Tigers/PCA/plinkFiles", pattern="*.gds", full.names=TRUE, recursive=FALSE)
# 
#for (f in files) {
#   print(f)
#   
#   gds = snpgdsOpen(f) #open genofile
#   sampIds = read.gdsn(index.gdsn(gds, "sample.id")) #grab sample ids 
#   famIds = as.data.frame(sampIds) %>%
#     left_join(popsDF, by = c("sampIds"="Individual")) %>%
#     pull(Subspecies2)#make family ids
#   
#   #LD prune because data set is small, set r2 = 0.7, maf 5%, and 50 snp window 
#   set.seed(10192021)
#   snpset = snpgdsLDpruning(gds, sample.id=sampIds, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F) 
#   pruned = unlist(snpset, use.names=FALSE)
#   
#   #Run KING
#   KING = snpgdsIBDKING(gds, sample.id = sampIds, snp.id=pruned, autosome.only = F)
#   
#   #Make GRM
#   KINGmat = KING$kinship %>%
#     as.matrix()  
#   colnames(KINGmat) = KING$sample.id 
#   rownames(KINGmat) = KING$sample.id 
#   snpgdsClose(gds)
#   
#   #Reopen GDS data 
#   genoFile = GdsGenotypeReader(filename = f)#read in GDS data
#   genoData = GenotypeData(genoFile)#create a GenotypeData class object
#   
#   ##Partition data into relateds and unrelated (less than first cousins)
#   sampset = pcairPartition(kinobj = KINGmat, kin.thresh=2^(-9/2), divobj = KINGmat, div.thresh=2^(-9/2))
#   unrelsID = sampset$unrels
#   
#   #####Look at relatedness in all
#   fullDataKinship = reshape2::melt(KINGmat) %>%
#     mutate(pop1 = popsDF$Subspecies2[match(Var1, popsDF$Individual)],
#            pop2 = popsDF$Subspecies2[match(Var2, popsDF$Individual)]) %>%
#     filter(pop1 == pop2 & Var1 != Var2)
#   
#   #####Look at relatedness in unrelated
#   unrelateds = reshape2::melt(KINGmat) %>%
#     filter(Var1%in%unrelsID & Var2%in%unrelsID) %>%
#     mutate(pop1 = popsDF$Subspecies2[match(Var1, popsDF$Individual)],
#            pop2 = popsDF$Subspecies2[match(Var2, popsDF$Individual)]) %>%
#     filter(pop1 == pop2 & Var1 != Var2)
# 
#   ####Make unrelateds files
#   unrelateds = as.data.frame(unrelsID) %>%
#     left_join(popsDF, by = c("unrelsID"="Individual")) %>%
#     select(unrelsID, Subspecies2)#make data frame
#   
#   fnameKin = gsub("PCA/plinkFiles", "Relatedness/FindDups/SNPRelate", f)
#   fnameUnrels = gsub("PCA/plinkFiles", "Relatedness/FindDups/SNPRelate/unrelateds_pcair", f)
# 
#   write.table(fullDataKinship, paste0(fnameKin,"_snpgdsIBDKING_Kinship.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
#   write.table(unrelateds, paste0(fnameUnrels,"_unrelateds.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# }

#####Compare kinship and truffle to find duplicates
files <- list.files(path="~/Documents/Tigers/Relatedness/FindDups/SNPRelate", pattern="*.txt", full.names=TRUE, recursive=FALSE)
truffleFiles_HC <- list.files(path="~/Documents/Tigers/Relatedness/FindDups/highcov-truffle-vcftools", pattern="*.refrmt", full.names=TRUE, recursive=FALSE)
truffleFiles_LC <- list.files(path="~/Documents/Tigers/Relatedness/FindDups/all-truffle-vcftools", pattern="*.refrmt", full.names=TRUE, recursive=FALSE)

#Keep only duplicates
df = rbindlist(sapply(files, read.delim, simplify = FALSE), use.names = TRUE, idcol = "FileName") %>% 
  filter(value > 0.49) %>%
  mutate(FileName = gsub(""))
  rename(ID1 = "Var1", ID2 = "Var2")

truf_highCov = rbindlist(sapply(truffleFiles_HC, read.delim, simplify = FALSE), use.names = TRUE, idcol = "FNameTruff") %>% 
  filter(IBD2 > 0.75) %>%
  left_join(df)
  
truf_lowCov = rbindlist(sapply(truffleFiles_LC, read.delim, simplify = FALSE), use.names = TRUE, idcol = "FNameTruff") %>% 
  filter(IBD2 > 0.75) %>%
  left_join(df)
