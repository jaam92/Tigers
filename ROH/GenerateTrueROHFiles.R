#load libraries
library(tidyverse)
library(data.table)

#Read ROH Files in 
setwd("~/Documents/Tigers/ROH")

df = read_delim("allChroms_allIndivs_garlicROH_700bpWindow_CoverageCallableSites.bed", delim = "\t", col_names = c("CHROM", "AUTO_START", "AUTO_END", "TYPE", "AUTO_LEN", "INDV",  "PropCovered"), col_types = "cnncncn")

##Identify ROH greater than or equal to 100Kb and remove ROH where the avg. prop covered by SNPs is within 1 SD of the mean 
ROHgr100kb = df[which(df$AUTO_LEN >= 100000),][c("TYPE","AUTO_START", "AUTO_END", "AUTO_LEN", "INDV", "PropCovered", "CHROM")]
z = data.table(ROHgr100kb)
z[,ToKeep := abs(ROHgr100kb$PropCovered - mean(ROHgr100kb$PropCovered)) < sd(ROHgr100kb$PropCovered)][ToKeep  == TRUE] #ID ROH to keep
FinalDF_allROH = subset(z, z$ToKeep == "TRUE") #Subset out true ROH

#Write final ROH to outfile
write.table(x = FinalDF_allROH, file = "TrueROH_propCoveredwithin1SDMean_gr100kb_allChroms_highCov_runSpeciesSep_garlic.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
