###This script will pull information from VCF file for each individual and pair it with corresponding VEP annotation at each site, then output a new file###
###Input files: VEP annotation file and Polarized VCF with ANCHOM,DERHOM,HET, or MISSING as GT
###Output file: File for each individual and chrom with CHROM, POS, GT, ANNOT, IMPACT, SIFT

#Load Libraries
library(data.table)
library(tidyverse)

#read in chroms
chrom = read.table("~/Documents/Tigers/AnnotSites/chromNames.txt") %>%
  pull(V1)
header = scan("~/Documents/Tigers/AnnotSites/VCFheader.txt", character(), quote = "")

for (k in chrom) {
  infile = paste0("~/Documents/Tigers/AnnotSites/PolarizedVCFs/AssignGT_CactusAncasREF_chr", k ,"_reformatted_masterFile_18SFS-ba-AN-MM-pcc-GM.vcf")
  skipHash = sprintf("grep -v '^#' %s", infile) #run bash command to remove hashes and take a variable for input file
  annotation = read_delim(file=paste0("~/Documents/Tigers/AnnotSites/Annotations/", k ,"_VEPandSIFT_felCat8_withTigerREFALT.txt"), delim="\t")
  vcf = fread(cmd=skipHash) %>% 
    select(-c(V1)) #remove the null column that get's added
  
  #Change column names
  colnames(vcf) = header
  names(vcf)[1] = "CHROM"
  names(vcf)[2] = "POS"
  colnames(vcf) = gsub(".*]","", colnames(vcf))
  
  #Add annotations
  annotatedVCF = vcf %>%
    filter(as.numeric(POS) %in% as.numeric(annotation$tigerPos)) %>%
    mutate(ANNOT = annotation$recodeConsequence[match(POS,annotation$tigerPos)],
           IMPACT = annotation$IMPACT[match(POS,annotation$tigerPos)],
           SIFT = annotation$SIFT_Consequence[match(POS,annotation$tigerPos)])
  rm(vcf) #make some space
  
  #Grab individuals
  indivs = colnames(annotatedVCF)[str_detect(colnames(annotatedVCF),"GT")]
 
  for (i in indivs){
    indivAnnotDF = annotatedVCF %>%
      select(CHROM, POS, i, ANNOT, IMPACT, SIFT)
    indivID = gsub(":GT", "", i)
    write.table(indivAnnotDF, file=paste0("~/Documents/Tigers/AnnotSites/AnnotatedVCF/", k, "_annotatedGTwithVEP_",indivID, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  cat(sprintf("done %s", k))
}
