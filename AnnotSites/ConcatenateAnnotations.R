#This script will concatenate all the annotation files either by individual or chromosome
library(tidyverse)
library(data.table)

#Load files
setwd("/scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/AnnotatedVCF")
fnames = list.files(pattern = "\\_addROHAnnot.txt$")

#Generate data frame
##create columns with fileName, population, and compute pi 
df = rbindlist(sapply(fnames, read.delim, col.names=c("CHROM", "POS", "GT", "ANNOT", "IMPACT", "SIFT", "withinTypeA", "withinTypeB", "withinTypeC"), simplify = FALSE), use.names = TRUE, idcol = "FileName") %>%
  mutate(ID = str_split_fixed(FileName, "_", 3)[,3],
         ID = gsub("_addROHAnnot.txt" ,"", ID),
         FileName = NULL)

#Make sure we have an annotation for every individual
individual_ids = read_delim("/scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/IndivFiles/individual_ids_unimputed_depthgrEqlTo5.txt", delim="\t")
table(individual_ids$Sample %in% unique(df$ID)) #should all be true

#split annotations by chromosome and individual
splitChroms = split(df, df$CHROM)
lapply(names(splitChroms), function(x){
  write.table(splitChroms[[x]], file = paste0("/scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/AnnotatedVCF/ConcatenatedFiles/allIndivs_", x, "_annotatedGTwithVEP.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)})

splitIndivs = split(df, df$ID)
lapply(names(splitIndivs), function(x){write.table(splitIndivs[[x]], file = paste0("/scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/AnnotatedVCF/ConcatenatedFiles/annotatedGTwithVEP_", x, "_allChroms.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)})
