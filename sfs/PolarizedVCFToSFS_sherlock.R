#Load Libraries
library(data.table)
library(tidyverse)

#Load files
setwd("/scratch/users/elliea/jazlyn-ellie/captive-tigers/sfs/Polarize")

fnames = list.files(pattern = "*_reformatted_masterFile_18SFS-ba-AN-MM-pcc-GM.vcf$")
pops = c("Generic","Amur","Bengal","Malayan","Sumatran")

for (pop in pops){
  indivs = read.delim("/scratch/users/elliea/jazlyn-ellie/captive-tigers/high-cov/18sAnd6s_annotated.txt") %>%
    filter(SubSpecies == pop)
  
  #Empty data frame to fill with summary data
  summaryInfo = data.frame()
  
  for (i in fnames){
    vcf = fread(file=fnames[i]) %>% 
      select(-c(V1)) #remove the null column that get's added
    chrom = str_split_fixed(fnames[i], "_", 4)[3] #grab chromosome 
    
    #Change column names
    colnames(vcf) = gsub(".*]","", colnames(vcf))
    colnames(vcf) = gsub(":GT", "", colnames(vcf))
    N18_indivs = colnames(vcf)[which(colnames(vcf) %in% indivs$Sample)]
    
    #Reformat and replace patterns with allele counts value
    GTCounts = vcf %>%
      select(all_of(N18_indivs)) %>%
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

  write.table(summaryInfo, file = paste(i, "_SFS_allChroms_SummaryFile_N18.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


}

