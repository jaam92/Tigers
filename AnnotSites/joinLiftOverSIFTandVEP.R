library(data.table)

#read in the vep annotation
vep = fread("grep -v '^##' ~/Documents/F2_variant.annot") %>%
  mutate(fel8pos = as.numeric(gsub(".*:","", Location)))

#read in positions from liftover
liftover = read_delim("~/Documents/fel5to8.txt", col_names = c("chrom", "fel5pos", "fel8pos"), delim = "\t")

#read in sift data then: 1) join by felcat5 pos to liftover; 2) join by felcat8 pos to vep
sift = fread("~/Documents/F2.gz") %>%
  rename(fel5pos=`#Position`) %>%
  left_join(liftover) %>% #attaches the chrom and felcat8 post
  left_join(vep) %>% #attaches all VEP annotations
  na.omit(cols="Consequence") #remove sift sites that don't have a vep annotation


finalDF = sift %>%
  select(chrom,fel8pos,Allele,SIFT_score,Consequence,IMPACT,Amino_acids,STRAND,Transcript_id,Gene_id,SYMBOL_SOURCE, SYMBOL,Feature_type) %>% #keep only cols of interest