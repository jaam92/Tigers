#Load Libraries
library(data.table)
library(tidyverse)

#Load files
setwd("~/Documents/TigerProject/Rscripts/sfs")
fnames = list.files(path="~/Documents/Tigers/AnnotSites/PolarizedVCFs", pattern="*pcc.vcf.gz", full.names=TRUE, recursive=FALSE) #input vcfs
pops = c("Generic","Amur","Bengal","Malayan","Sumatran")

for (pop in pops){
  indivs = read.delim("~/Documents/Tigers/IndivFiles/N10AndN6_unrelateds.txt") %>%
    filter(Subspecies == pop & N10 == 1)
  
  #Empty data frame to fill with summary data
  summaryInfo = data.frame()
  
  for (i in seq_along(fnames)){
    chrom = str_split_fixed(fnames[i], "_", 4)[3] #grab chromosome
    vcf = fread(file= fnames[i], sep = "\t",  fill = T) #read in vcf. fread does not like the hash and it forces you to add an extra column
    
    #Modify column names
    newColNames = colnames(vcf)[-length(colnames(vcf))] #grab column names and remove the extra column that gets added 
    newColNames = gsub(".*]","", newColNames)
    newColNames = gsub(":GT", "", newColNames)
    
    #After making column names re-add names and remove extra column from reading vcf in
    vcf = vcf %>% 
      select(-c("#[1]CHROM")) #remove the null column that gets added
    colnames(vcf) = newColNames
    sfs_indivs = colnames(vcf)[which(colnames(vcf) %in% indivs$ID)]
    
    #Reformat and replace patterns with allele counts value
    GTCounts = vcf %>%
      select(all_of(sfs_indivs)) %>%
      filter_all(all_vars(. != "Missing")) %>% #remove any sites with missing data
      group_by_all() %>% 
      count() %>% #count of all patterns in rows
      ungroup() %>% 
      mutate(across(everything(), ~replace(., . ==  "DerHom" , 2))) %>% 
      mutate(across(everything(), ~replace(., . ==  "Het" , 1))) %>% 
      mutate(across(everything(), ~replace(., . ==  "AncHom" , 0))) %>%
      rename("Count" = n)
    
    rm(vcf)
    
    #count of the frequency of each pattern
    FreqBin = GTCounts %>% 
      select(-c(Count)) %>% 
      mutate_all(~parse_number(.x)) %>% 
      rowSums() 
    
    #Data frame with frequency and counts
    Polarized = cbind.data.frame(FreqBin, GTCounts$Count) %>% 
      group_by(FreqBin) %>% 
      summarise(sum = sum(`GTCounts$Count`)) %>%
      mutate(chromosome = chrom)
      
    #Data frame with all chromosomes
    summaryInfo = bind_rows(summaryInfo, Polarized)
  
  }

  write.table(summaryInfo, file = paste0(pop, "_SFS_allChroms_SummaryFile_N10.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


}