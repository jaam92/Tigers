####This script will take in annotation data across all chromosome for each individual and count up the types of annotations (with three different counting methods), how many annotations fall within or outside of ROH of a given length, compute OR of NS to SYN variants within vs outside of ROH for each individual and across populations

####Input file: annotation file for each individual with data from all chroms
####Output file: File with above information for all individuals and a file with OR for different counting methods split by population 

#Load Libraries
library(data.table)
library(tidyverse)

#function to compute counts
summarizeCounts = function(dataFrame, col_name){
  colOfInterest = enquo(col_name)
  df = dataFrame %>%
    filter(GT != "AncHom" & GT != "Missing") %>% #Remove AncHom and Missing
    group_by(GT,ANNOT, !!colOfInterest) %>%
    count(name = "CountVariants") %>% #Count Variants (DERHOM =1 , HET = 1)
    mutate(CountAlleles = ifelse(GT == "DerHom", as.numeric(CountVariants)*2, as.numeric(CountVariants)),#Count Alleles (DERHOM =2 , HET = 1)
           CountDerHom = ifelse(GT == "DerHom", as.numeric(CountVariants), as.numeric(0))) %>% #Count DerHOM (DERHOM =1 )
    ungroup() %>%
    group_by(ANNOT, !!colOfInterest) %>%
    summarise_at(c("CountVariants", "CountAlleles", "CountDerHom"), sum)  %>%
    ungroup() 
  return(df)
}

summarizeCounts_SIFT = function(dataFrame, col_name){
  colOfInterest = enquo(col_name)
  df = dataFrame %>%
    na.omit(SIFT) %>%
    filter(GT != "AncHom" & GT != "Missing") %>% #Remove AncHom and Missing
    mutate(SIFTAnnot = ifelse(ANNOT == "NS" & SIFT == "deleterious", "PutDel", "PutNeu")) %>% 
    group_by(GT, SIFTAnnot, !!colOfInterest) %>%
    count(name = "CountVariants") %>% #Count Variants (DERHOM =1 , HET = 1)
    mutate(CountAlleles = ifelse(GT == "DerHom", as.numeric(CountVariants)*2, as.numeric(CountVariants)),#Count Alleles (DERHOM =2 , HET = 1)
           CountDerHom = ifelse(GT == "DerHom", as.numeric(CountVariants), as.numeric(0))) %>% #Count DerHOM (DERHOM =1 )
    ungroup() %>%
    group_by(SIFTAnnot, !!colOfInterest) %>%
    summarise_at(c("CountVariants", "CountAlleles", "CountDerHom"), sum)  %>%
    ungroup() 
  return(df)
}

pivotSummarizedCounts = function(dataFrame, col_name, TypeROH){
  colOfInterest = enquo(col_name)
  df = dataFrame %>% 
    pivot_longer(cols = starts_with("Count"), names_to = "Method") %>%
    arrange(ANNOT) %>%
    mutate(col_name  = ifelse(!!colOfInterest == 1, "ROH", "nonROH"),
           NewColName = paste(ANNOT,Method,col_name,TypeROH, sep = "_")) %>%
    select(NewColName,value) %>%
    spread(NewColName,value)
  return(df)
}

pivotSummarizedCounts_SIFT = function(dataFrame, col_name, TypeROH){
  colOfInterest = enquo(col_name)
  df = dataFrame %>% 
    pivot_longer(cols = starts_with("Count"), names_to = "Method") %>%
    arrange(SIFTAnnot) %>%
    mutate(col_name  = ifelse(!!colOfInterest == 1, "ROH", "nonROH"),
           NewColName = paste(SIFTAnnot,Method,col_name,TypeROH, sep = "_")) %>%
    select(NewColName,value) %>%
    spread(NewColName,value)
  return(df)
}

#Get input files
setwd("/scratch/users/elliea/jazlyn-ellie/captive-tigers/final_files/AnnotsVEPandSIFT/annotPolarizedVCF/annot_ROHVEPSIFT/AllChroms")
#setwd("~/Documents/Tigers/AnnotSites/AnnotatedVCF/AllChroms/")#run locally

individual_ids = read_delim("/scratch/users/elliea/jazlyn-ellie/captive-tigers/final_files/SampleLists/TableX-SampleDetails.txt", delim = "\t") %>%
  filter(`Depth (post filtering)` >= 5 & Coverage_group == "Unimputed" & is.na(Removed))#meta-data

#individual_ids = read_delim("~/Documents/Tigers/IndivFiles/TableX-SampleDetails.txt", delim = "\t") %>%
#  filter(`Depth (post filtering)` >= 5 & Coverage_group == "Unimputed" & is.na(Removed))#meta-data run locally

filenames = list.files(pattern = glob2rx("annotatedGTwithVEP_*_allChroms.txt"))

#Empty data frame to fill with data
allIndiv = list()

#Loop thru all individuals files
for (i in 1:length(filenames)){
  
  #Read files in 
  annots = read.delim(file = filenames[i])
  indiv = str_match(filenames[i], "annotatedGTwithVEP_\\s*(.*?)\\s*_allChroms.txt")[,2] #pull individual ID from between underscores
  

####Grab counts ROH for Type A
  CountDF_TypeA = summarizeCounts(annots, TypeA)
  CountDF_TypeA_SIFT = summarizeCounts_SIFT(annots, TypeA)

####Grab counts ROH for Type B
  CountDF_TypeB = summarizeCounts(annots, TypeB)
  CountDF_TypeB_SIFT = summarizeCounts_SIFT(annots, TypeB)
  
####Grab counts ROH for Type C
  CountDF_TypeC = summarizeCounts(annots, TypeC)
  CountDF_TypeC_SIFT = summarizeCounts_SIFT(annots, TypeC)
  
####Grab counts within and outside Type A ROH
  ROHStats_TypeA = pivotSummarizedCounts(CountDF_TypeA, TypeA, "TypeA")
  ROHStats_TypeA_SIFT = pivotSummarizedCounts_SIFT(CountDF_TypeA_SIFT, TypeA, "TypeA")
    
####Grab counts within and outside Type B ROH
  ROHStats_TypeB = pivotSummarizedCounts(CountDF_TypeB, TypeB, "TypeB")
  ROHStats_TypeB_SIFT = pivotSummarizedCounts_SIFT(CountDF_TypeB_SIFT, TypeB, "TypeB")
  
####Grab counts within and outside Type C ROH
  ROHStats_TypeC = pivotSummarizedCounts(CountDF_TypeC, TypeC, "TypeC")
  ROHStats_TypeC_SIFT = pivotSummarizedCounts_SIFT(CountDF_TypeC_SIFT, TypeC, "TypeC")

###Counts of just annotations it will be the same across any df
  CountDF_withJustAnnots = CountDF_TypeA %>% 
    group_by(ANNOT) %>% 
    summarise_at(.vars= vars(CountVariants, CountAlleles, CountDerHom), .funs = sum) 
  
####Counts of variants by Impact
  ImpactCounts = annots %>%
    group_by(IMPACT) %>%
    count() %>%
    spread(IMPACT,n)
  
####Count up the total number of annotated variants and how they are distributed
  Calls = annots %>%
    group_by(GT) %>%
    count() %>%
    spread(GT,n) %>%
    mutate(LineCount = dim(annots)[1])

####Counts using SIFT
  SIFTAnnot = annots %>% 
    na.omit(SIFT) %>%
    filter(GT != "AncHom", GT != "Missing") %>%
    mutate(SIFTAnnot = ifelse(ANNOT == "NS" & SIFT == "deleterious", "PutDel", "PutNeu")) %>% 
    group_by(GT,SIFTAnnot) %>%
    count(name = "CountVariants") %>%
    mutate(CountAlleles = ifelse(GT == "DerHom", as.numeric(CountVariants)*2, as.numeric(CountVariants)),
           CountDerHom = ifelse(GT == "DerHom", as.numeric(CountVariants), as.numeric(0))) %>%
    group_by(SIFTAnnot) %>%
    summarise_at(c("CountVariants", "CountAlleles", "CountDerHom"), sum) %>%
    ungroup() %>%
    pivot_longer(cols = starts_with("Count"), names_to = "Method") %>% #convert table from wide to long
    mutate(NewColName = paste0(SIFTAnnot,"SIFT_",Method)) %>%
    select(NewColName,value) %>% #only keep new col and values
    spread(NewColName,value)

####Put all the counts together and add odds-ratios within vs outside ROH  
  indivWideDF = CountDF_withJustAnnots %>% #convert table from wide to long
    pivot_longer(cols = starts_with("Count"), names_to = "Method") %>%
    mutate(NewColName = paste0(ANNOT,"_",Method)) %>%
    select(NewColName,value) %>% #only keep new col and values
    spread(NewColName,value) %>% #make data one row
    cbind.data.frame(Calls, ImpactCounts, SIFTAnnot, ROHStats_TypeA, ROHStats_TypeB, ROHStats_TypeC, ROHStats_TypeA_SIFT, ROHStats_TypeB_SIFT, ROHStats_TypeC_SIFT) %>% #add on total calls and impact counts
    mutate(ID = indiv, 
           Subspecies = individual_ids$Subspecies_GroupID_Corrected[match(ID, individual_ids$Sample)],
           OR_AlleleCopies_TypeA = tryCatch({(SY_CountAlleles_nonROH_TypeA*NS_CountAlleles_ROH_TypeA)/(SY_CountAlleles_ROH_TypeA*NS_CountAlleles_nonROH_TypeA)}, error=function(e) NA),
           OR_Variants_ROH_TypeA = tryCatch({(SY_CountVariants_nonROH_TypeA*NS_CountVariants_ROH_TypeA)/(SY_CountVariants_ROH_TypeA*NS_CountVariants_nonROH_TypeA)}, error=function(e) NA),
           OR_DerHom_ROH_TypeA = tryCatch({(SY_CountDerHom_nonROH_TypeA*NS_CountDerHom_ROH_TypeA)/(SY_CountDerHom_ROH_TypeA*NS_CountDerHom_nonROH_TypeA)}, error=function(e) NA),
           OR_AlleleCopies_TypeB = tryCatch({(SY_CountAlleles_nonROH_TypeB*NS_CountAlleles_ROH_TypeB)/(SY_CountAlleles_ROH_TypeB*NS_CountAlleles_nonROH_TypeB)}, error=function(e) NA),
           OR_Variants_ROH_TypeB = tryCatch({(SY_CountVariants_nonROH_TypeB*NS_CountVariants_ROH_TypeB)/(SY_CountVariants_ROH_TypeB*NS_CountVariants_nonROH_TypeB)}, error=function(e) NA),
           OR_DerHom_ROH_TypeB = tryCatch({(SY_CountDerHom_nonROH_TypeB*NS_CountDerHom_ROH_TypeB)/(SY_CountDerHom_ROH_TypeB*NS_CountDerHom_nonROH_TypeB)}, error=function(e) NA),
           OR_AlleleCopies_TypeC = tryCatch({(SY_CountAlleles_nonROH_TypeC*NS_CountAlleles_ROH_TypeC)/(SY_CountAlleles_ROH_TypeC*NS_CountAlleles_nonROH_TypeC)}, error=function(e) NA),
           OR_Variants_ROH_TypeC = tryCatch({(SY_CountVariants_nonROH_TypeC*NS_CountVariants_ROH_TypeC)/(SY_CountVariants_ROH_TypeC*NS_CountVariants_nonROH_TypeC)}, error=function(e) NA),
           OR_DerHom_ROH_TypeC = tryCatch({(SY_CountDerHom_nonROH_TypeC*NS_CountDerHom_ROH_TypeC)/(SY_CountDerHom_ROH_TypeC*NS_CountDerHom_nonROH_TypeC)}, error=function(e) NA)) %>% #add Individual ID, Population, and OR for counting methods and different types of ROH
    select(ID, Subspecies, everything())#move individual id and population to front
  
  #Dataframe with all individuals (use bind rows so columns w/out values filled as na)
  allIndiv[[i]] = indivWideDF
  
}

allIndivDF = bind_rows(allIndiv) %>%
          mutate_if(is.numeric, round, digits=4)

#Write output file
write.table(allIndivDF, file="/scratch/users/elliea/jazlyn-ellie/captive-tigers/final_files/AnnotsVEPandSIFT/annotPolarizedVCF/GTAnnotationCountResults_Nov2022_Tigers_addSIFT.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#write.table(allIndivDF, file="~/Documents/Tigers/AnnotSites/GTAnnotationCountResults_Nov2022_Tigers_addSIFT.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE) #run locally

rm(allIndiv)

#Calculate OR for entire population
PerPopulationOR = allIndivDF %>%
  group_by(Subspecies) %>%
  summarise(OR_AlleleCopies_TypeA =(sum(SY_CountAlleles_nonROH_TypeA, na.rm = TRUE)*sum(NS_CountAlleles_ROH_TypeA, na.rm = TRUE))/(sum(SY_CountAlleles_ROH_TypeA, na.rm = TRUE)*sum(NS_CountAlleles_nonROH_TypeA, na.rm = TRUE)),
            OR_Variants_TypeA = (sum(SY_CountVariants_nonROH_TypeA, na.rm = TRUE)*sum(NS_CountVariants_ROH_TypeA, na.rm = TRUE))/(sum(SY_CountVariants_ROH_TypeA, na.rm = TRUE)*sum(NS_CountVariants_nonROH_TypeA, na.rm = TRUE)),
            OR_DerHom_TypeA = (sum(SY_CountDerHom_nonROH_TypeA, na.rm = TRUE)*sum(NS_CountDerHom_ROH_TypeA, na.rm = TRUE))/(sum(SY_CountDerHom_ROH_TypeA, na.rm = TRUE)*sum(NS_CountDerHom_nonROH_TypeA, na.rm = TRUE)),
            OR_AlleleCopies_TypeB =(sum(SY_CountAlleles_nonROH_TypeB, na.rm = TRUE)*sum(NS_CountAlleles_ROH_TypeB, na.rm = TRUE))/(sum(SY_CountAlleles_ROH_TypeB, na.rm = TRUE)*sum(NS_CountAlleles_nonROH_TypeB, na.rm = TRUE)),
            OR_Variants_TypeB = (sum(SY_CountVariants_nonROH_TypeB, na.rm = TRUE)*sum(NS_CountVariants_ROH_TypeB, na.rm = TRUE))/(sum(SY_CountVariants_ROH_TypeB, na.rm = TRUE)*sum(NS_CountVariants_nonROH_TypeB, na.rm = TRUE)),
            OR_DerHom_TypeB = (sum(SY_CountDerHom_nonROH_TypeB, na.rm = TRUE)*sum(NS_CountDerHom_ROH_TypeB, na.rm = TRUE))/(sum(SY_CountDerHom_ROH_TypeB, na.rm = TRUE)*sum(NS_CountDerHom_nonROH_TypeB, na.rm = TRUE)),
            OR_AlleleCopies_TypeC =(sum(SY_CountAlleles_nonROH_TypeC, na.rm = TRUE)*sum(NS_CountAlleles_ROH_TypeC, na.rm = TRUE))/(sum(SY_CountAlleles_ROH_TypeC, na.rm = TRUE)*sum(NS_CountAlleles_nonROH_TypeC, na.rm = TRUE)),
            OR_Variants_TypeC = (sum(SY_CountVariants_nonROH_TypeC, na.rm = TRUE)*sum(NS_CountVariants_ROH_TypeC, na.rm = TRUE))/(sum(SY_CountVariants_ROH_TypeC, na.rm = TRUE)*sum(NS_CountVariants_nonROH_TypeC, na.rm = TRUE)),
            OR_DerHom_TypeC = (sum(SY_CountDerHom_nonROH_TypeC, na.rm = TRUE)*sum(NS_CountDerHom_ROH_TypeC, na.rm = TRUE))/(sum(SY_CountDerHom_ROH_TypeC, na.rm = TRUE)*sum(NS_CountDerHom_nonROH_TypeC, na.rm = TRUE))) %>%
          mutate_if(is.numeric, round, digits=4) 

#Write output file
#write.table(PerPopulationOR, file="/scratch/users/elliea/jazlyn-ellie/captive-tigers/final_files/AnnotsVEPandSIFT/annotPolarizedVCF/PerPopulationOR_Nov2022_Tigers.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#write.table(PerPopulationOR, file="~/Documents/Tigers/AnnotSites/PerPopulationOR_Nov2022_Tigers.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE) #run locally



