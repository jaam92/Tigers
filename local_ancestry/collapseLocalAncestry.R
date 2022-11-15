library(data.table)
library(GenomicRanges)
library(tidyverse)


#set working directory and read files in
setwd("~/Tigers/local_ancestry/local_calls/")
fnames = list.files(pattern = "\\.msp.tsv")

#Function to collapse local ancestry 
collapseRanges = function(ancestryDataFrame, ancestryGroup){
  df = local_anc %>% 
    filter(ancestry == ancestryGroup) 
  gr = GRanges(df) #https://www.biostars.org/p/386585/#423870
  final = unlist(reduce(split(gr, gr$sample))) #https://www.biostars.org/p/386585/#423870
  sample = names(final)
  final_df = final %>%
    as_tibble() %>%
    mutate(sampleID = sample,
           ancestry = ancestryGroup,
           strand = NULL) %>%
    dplyr::rename(length = width) 
  return(final_df)
}


#Read in local ancestry and rename columns
local_anc = rbindlist(sapply(fnames, fread, simplify = FALSE), use.names = TRUE) %>%
  select(-c(sgpos,egpos,`n snps`))  %>%
  pivot_longer(cols = 4:141, names_to = "sample", values_to = "ancestry")

colnames(local_anc) = c("chrom", "start", "end", "sample","ancestry")

#Collapse ranges per individual and ancestry from
datalist = list()
datalist = vector("list", length = 5) # or pre-allocate for slightly more efficiency

for (anc in 0:4) {
  tempDF = collapseRanges(local_anc, anc)
  i = anc + 1
  datalist[[i]] <- tempDF # add it to your list
}

all_data = do.call(rbind, datalist)
write.table(all_data, file = "allChroms_localAncestry_perIndiv.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")