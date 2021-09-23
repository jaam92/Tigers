library(data.table)

#read in the vep annotation
vep = fread("grep -v '^##' ~/Documents/Tigers/AnnotSites/F2_variant.annot") %>%
  mutate(fel8pos = as.numeric(gsub(".*:","", Location))) %>%
  filter(IMPACT != "MODIFIER") %>% #remove intergenic stuff
  mutate(IMPACT = str_replace_all(IMPACT, c("LOW" = "0", "MODERATE" = "1", "HIGH" = "2"))) %>% #make heirarchy numeric so we can select the most damaging annotation for a given transcript
  group_by(fel8pos) %>% 
  slice(which.max(as.numeric(IMPACT))) %>% #choose most damaging annot and if there is a tie for impact pick first reported (normally ends up being the canonical)
  mutate(IMPACT = str_replace_all(IMPACT, c("0" = "LOW", "1" = "MODERATE", "2" = "HIGH")), #change the impact column back to original annot
         REF_ALT = str_sub(`#Uploaded_variation`, start= -3)) %>%  
  ungroup() %>%
  mutate(recodeConsequence = str_replace_all(Consequence, c("missense_variant,splice_region_variant" = "NS", "splice_donor_variant,missense_variant" = "NS",  "splice_region_variant,synonymous_variant" = "SY", "splice_region_variant,synonymous_variant" = "SY", "stop_retained_variant" = "SY",  "stop_gained,splice_region_variant" = "LOF", "stop_lost" = "LOF", "start_lost" = "LOF", "start_lost,synonymous_variant" = "LOF", "stop_gained,start_lost" = "LOF", "missense_variant" = "NS", "stop_gained" = "LOF", "synonymous_variant" = "SY"))) %>% #recode to match script and categorize as either NS, SY, LOF
  filter(!str_detect(recodeConsequence, 'splice')) %>% #remove splice_acceptor, splice_donor, splice_region annotations
  separate_rows(Amino_acids, REF_ALT, sep="/")

         
#read in positions from liftover
liftover = read_delim("~/Documents/Tigers/AnnotSites/fel5to8.txt", col_names = c("chrom", "fel5pos", "fel8pos"), delim = "\t")

#read in sift data then: 1) join by felcat5 pos to liftover; 2) join by felcat8 pos to vep
sift = fread("~/Documents/Tigers/AnnotSites/F2.gz") %>%
  rename(fel5pos=`#Position`) %>%
  left_join(liftover)  #attaches the chrom and felcat8 post
  

finalDF = vep %>%
  left_join(sift, by = c("fel8pos" = "fel8pos", "Amino_acids"="New_amino_acid", "REF_ALT" = "New_allele")) %>% #attaches all VEP annotations
  na.omit(cols="Consequence") %>% #remove sift sites that don't have a vep annotation
  select(chrom,fel8pos,REF_ALT,SIFT_score,Consequence,recodeConsequence,IMPACT,Amino_acids,STRAND,Transcript_id,Gene_id,SYMBOL_SOURCE, SYMBOL,Feature_type) #keep only cols of interest
