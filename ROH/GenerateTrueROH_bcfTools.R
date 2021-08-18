#Load libraries and data
library(dplyr)
library(data.table)

#Load bcftools roh outputs and meta data
setwd("/u/scratch/j/jmooney3/bigCats/bcftoolsROH/Ranges")

fnames = list.files(pattern = "\\_rohs.txt$")

nonQCROH = rbindlist(sapply(fnames, fread, simplify = FALSE), use.names = TRUE, idcol = "FileName") %>%
  mutate(LengthROH = as.numeric(`[3]End`) - as.numeric(`[2]Start`),
         Sample = gsub("final_","",FileName),
         Sample = gsub("_rohs.txt","",Sample),
         FileName = NULL) %>%
  filter(LengthROH >= 100000) #minimum run length of 100Kb 

#QC all roh greater than 50kb
z = data.table(nonQCROH)
z[,ToKeep := abs(nonQCROH$`[5]MeanQual` - mean(nonQCROH$`[5]MeanQual`)) < sd(nonQCROH$`[5]MeanQual`)][ToKeep  == TRUE] #create variable that identifies what to drop or keep based on mean quality
FinalDF_allROH = subset(z, z$ToKeep == "TRUE")

#Write final ROH to outfile
write.table(x = FinalDF_allROH,file = "TrueROH_1SDMeanQual_gr100kb_allChroms_highCov_runSpeciesSep_bcfTools.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
