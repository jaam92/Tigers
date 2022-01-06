library(tidyverse)


for (i in 1:18) {
  
  posIndex = read_delim(file = paste0("/scratch/users/elliea/jazlyn-ellie/captive-tigers/addCallableSites_ROHandIBD/vcfSites/", i, "_posWithRowIndex.out"), delim = "\t", col_names = c("chrom", "pos", "index"), col_types = "cnn") 
  
  ibd = read_delim(file = paste0("/scratch/users/elliea/jazlyn-ellie/captive-tigers/addCallableSites_ROHandIBD/IBD/", i, "_truffle_allSubSpecies_calledPerSpecies.segments.bed"), delim = "\t", col_names = c("chrom","startVar","endVar", "sample1","sample2", "typeIBD", "IBDLengthMb"), col_types = "cnncccn") %>% 
    mutate(startPos = posIndex$pos[match(startVar,posIndex$index)],
           endPos = posIndex$pos[match(endVar,posIndex$index)],
           IBDLengthMbFromPos = (endPos-startPos)/10^6) 
  
  reformatBedFile = ibd %>%
    select("chrom", "startPos", "endPos", "sample1","sample2", "typeIBD", "IBDLengthMb")
  
  write.table(reformatBedFile, paste0("/scratch/users/elliea/jazlyn-ellie/captive-tigers/addCallableSites_ROHandIBD/IBD/", i, "_truffle_allSubSpecies_calledPerSpecies.segments.convert2VCFpos.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  print(i)
}




