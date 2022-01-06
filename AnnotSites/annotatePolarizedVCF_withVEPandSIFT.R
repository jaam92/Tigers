###This script will pull information from VCF file for each individual and pair it with corresponding VEP annotation at each site, then output a new file###
###Input files: VEP annotation file and Polarized VCF with ANCHOM,DERHOM,HET, or MISSING as GT
###Output file: File for each individual and chrom with CHROM, POS, GT, ANNOT, IMPACT, SIFT

#Load Libraries
library(data.table)
library(tidyverse)

#read in chroms
chrom = read.table("~/Documents/Tigers/AnnotSites/chromNames.txt") %>%
  pull(V1)

for (j in chrom) {
  #Load file
  annotation = read_delim(file=paste0("~/Documents/Tigers/AnnotSites/Annotations/", j ,"_VEPandSIFT_felCat8_withTigerREFALT.txt"), delim="\t")
  
  infileVCF = fread(file = paste0("~/Documents/Tigers/AnnotSites/PolarizedVCFs/AssignGT_CactusAncasREF_chr", j ,"_reformatted_masterFile_18SFS-ba-AN-MM-pcc-GM.vcf"), sep = "\t",  fill = T) #read in vcf fread does not like the hash and it forces you to add an extra column
  
  #Modify column names
  newColNames = colnames(infileVCF)[-length(colnames(infileVCF))] #grab column names and remove the extra column that gets added 
  newColNames = gsub(".*]","", newColNames)
  
  #After making column names re-add names and remove extra column from reading vcf in
  infileVCF = infileVCF %>% 
    select(-c("#[1]CHROM")) #remove the null column that gets added
  colnames(infileVCF) = newColNames
  
  #Add annotations
  annotatedVCF = infileVCF  %>%
    filter(as.numeric(POS) %in% as.numeric(annotation$tigerPos)) %>%
    mutate(ANNOT = annotation$recodeConsequence[match(POS,annotation$tigerPos)],
           IMPACT = annotation$IMPACT[match(POS,annotation$tigerPos)],
           SIFT = annotation$SIFT_Consequence[match(POS,annotation$tigerPos)])
  rm(vcf) #make some space
  
  #Grab individuals
  indivs = colnames(annotatedVCF)[str_detect(colnames(annotatedVCF),"GT")]
 
  for (i in indivs){
    indivAnnotDF = annotatedVCF %>%
      select(CHROM, POS, all_of(i), ANNOT, IMPACT, SIFT)
    colnames(indivAnnotDF) = gsub(":GT", "", colnames(indivAnnotDF))
    indivID = gsub(":GT", "", i)
    write.table(indivAnnotDF, file=paste0("~/Documents/Tigers/AnnotSites/AnnotatedVCF/", j, "_annotatedGTwithVEP_",indivID, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  cat(sprintf("done %s", j))
}
