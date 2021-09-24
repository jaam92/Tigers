#load libraries
library(tidyverse)
library(data.table)
library(optparse)

###Use input parameters from command line, if not found set to default
option_list <- list( 
  make_option(c("--VEP"), type="character", default="/scratch/users/elliea/jazlyn-ellie/captive-tigers/VEP/annotsVEP/F2_variant.annot", 
              help="annotations from VEP for tiger using felcat8 coords", metavar="character"),
  make_option(c("--SIFT"), type="character", default="/scratch/users/elliea/jazlyn-ellie/captive-tigers/liftover/Felis_catus_6.2.83/F2.gz", 
              help="felcat5 SIFT annotations from ensembl", metavar="character"),
  make_option(c("--LiftOverCoords"), type="character", default="/scratch/users/elliea/jazlyn-ellie/captive-tigers/liftover/runJazlynScript/F2_LiftOverfelCat8_combo.txt", 
              help="liftover coords from felcat5 sift to felcat8", metavar="character"),
  make_option(c("--outFile"), type="character", default="/scratch/users/elliea/jazlyn-ellie/captive-tigers/AnnotsVEPandSIFT/F2_VEPandSIFT_felCat8_withTigerREFALT.txt",
              help="specify output file name with full path", metavar="character")
)

#Parse the parameters and assign to variables
opt = parse_args(OptionParser(option_list=option_list)) #this must come first before we assign 
infileVEP = opt$VEP
infileLiftOver = opt$LiftOverCoords
infileSIFT = opt$SIFT
outputFile = opt$outFile

###Start the work

#read in the vep annotation
skipHash = sprintf("grep -v '^##' %s", infileVEP) #run bash command to remove hashes and take a variable for input file
#vep = fread("grep -v '^##' ~/Documents/Tigers/AnnotSites/F2_variant.annot") %>% #run local
vep = fread(cmd=skipHash) %>%
  mutate(fel8pos = as.numeric(gsub(".*:","", Location))) %>%
  filter(IMPACT != "MODIFIER") %>% #remove intergenic stuff
  mutate(IMPACT = str_replace_all(IMPACT, c("LOW" = "0", "MODERATE" = "1", "HIGH" = "2"))) %>% #make heirarchy numeric so we can select the most damaging annotation for a given transcript
  group_by(fel8pos) %>% 
  slice(which.max(as.numeric(IMPACT))) %>% #choose most damaging annot and if there is a tie for impact pick first reported (normally ends up being the canonical)
  mutate(IMPACT = str_replace_all(IMPACT, c("0" = "LOW", "1" = "MODERATE", "2" = "HIGH")), #change the impact column back to original annot
         REF_ALT = as.character("REF/ALT")) %>%  
  ungroup() %>%
  mutate(recodeConsequence = str_replace_all(Consequence, c("missense_variant,splice_region_variant" = "NS", "splice_donor_variant,missense_variant" = "NS",  "splice_region_variant,synonymous_variant" = "SY", "splice_region_variant,synonymous_variant" = "SY", "stop_retained_variant" = "SY",  "stop_gained,splice_region_variant" = "LOF", "stop_lost" = "LOF", "start_lost" = "LOF", "start_lost,synonymous_variant" = "LOF", "stop_gained,start_lost" = "LOF", "missense_variant" = "NS", "stop_gained" = "LOF", "synonymous_variant" = "SY"))) %>% #recode to match script and categorize as either NS, SY, LOF
  separate_rows(Amino_acids, REF_ALT, sep="/") %>%
  filter(!str_detect(recodeConsequence, 'splice') & REF_ALT == "ALT") #remove splice_acceptor, splice_donor, splice_region annotations

#read in positions from liftover
#liftover = read_delim("~/Documents/Tigers/AnnotSites/F2_LiftOverfelCat8_combo.txt", col_names = c("chrom", "fel5pos", "fel8pos"), delim = "\t") #run local
liftover = read_delim(infileLiftOver, col_names = c("chrom", "fel5pos", "fel8pos"), delim = "\t")

#read in sift data then: 1) join by felcat5 pos to liftover; 2) join by felcat8 pos to vep
#sift = fread("~/Documents/Tigers/AnnotSites/F2.gz") %>% #run local
sift = fread(infileSIFT) %>%
  rename(fel5pos=`#Position`) %>%
  left_join(liftover)  #attaches the chrom and felcat8 post

finalDF = vep %>%
  left_join(sift, by = c("fel8pos" = "fel8pos", "Amino_acids"="New_amino_acid")) %>% #attaches all VEP annotations
  na.omit(cols="Consequence") %>% #remove sift sites that don't have a vep annotation
  mutate(REF_ALT_alleles = str_sub(`#Uploaded_variation`, start= -3)) %>%
  mutate(SIFT_Consequence = ifelse(recodeConsequence != "SY" & as.numeric(SIFT_score) < 0.05, "deleterious", "benign")) %>%
  select(chrom,fel8pos,REF_ALT_alleles,Allele,SIFT_score,SIFT_Consequence,Consequence,recodeConsequence,IMPACT,Amino_acids,STRAND,Transcript_id,Gene_id,SYMBOL_SOURCE, SYMBOL,Feature_type) %>% #keep only cols of interest
  distinct(.keep_all = TRUE) #remove duplicate rows that occur when join with sift

#write output and clean up 
#write.table(finalDF, "~/Documents/Tigers/AnnotSites/VEPandSIFT_felCat8_withTigerREFALT.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE) #run local
write.table(finalDF, outputFile, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE) 
