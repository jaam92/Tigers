```
library(tidyverse)
library(dplyr)
library(reshape2)
library(ggpubr)
library(corrplot)
library(data.table)
library(tidyverse)
library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(corrplot)
library(ggdendro)
library(dendextend)
```
```
fnames = list.files(pattern = "\\.relatedness1")
gnames = list.files(pattern = "\\.relatedness2")
hnames = list.files(pattern = "\\_highcov.csv$")
```
#Generate data frames
##create columns with fileName, population, and compute pi 
```
relatedness_df = rbindlist(sapply(fnames, read.delim, simplify = FALSE), use.names = TRUE, idcol = "FileName")
relatedness2_df = rbindlist(sapply(gnames, read.delim, simplify = FALSE), use.names = TRUE, idcol = "FileName")
ibd_mle_df = rbindlist(sapply(hnames, read.csv, simplify = FALSE), use.names = TRUE, idcol = "FileName")
```
# clean files for plotting
```
relatedness_df <- relatedness_df %>%
  rowwise() %>%
  mutate(unique_id = paste(sort(c(INDV1, INDV2)), collapse = '-')) %>%
  ungroup() %>%
  subset(INDV1 != INDV2) %>%
  mutate(FileName = replace(FileName, FileName == "amurs-highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz.relatedness1.relatedness", "Amur")) %>%
  mutate(FileName = replace(FileName, FileName == "generic-highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz.relatedness1.relatedness", "Generic")) %>%
  mutate(FileName = replace(FileName, FileName == "bengals-highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz.relatedness1.relatedness", "Bengal")) %>%
  mutate(FileName = replace(FileName, FileName == "malayan-highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz.relatedness1.relatedness", "Malayan")) %>%
  mutate(FileName = replace(FileName, FileName == "sumatran-highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz.relatedness1.relatedness", "Sumatran")) %>%
  mutate(FileName = replace(FileName, FileName == "indochinese-highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz.relatedness1.relatedness", "Indochinese")) %>%
  select(-c(INDV1,INDV2)) %>%
  rename(Subspecies=FileName)

relatedness2_df <- relatedness2_df %>%
  rowwise() %>%
  mutate(unique_id = paste(sort(c(INDV1, INDV2)), collapse = '-')) %>%
  ungroup() %>%
  subset(INDV1 != INDV2) %>%
  select(-c(FileName,N_AaAa,N_AAaa,N1_Aa,N2_Aa)) 

#reorder and sample every other to get rid of duplicates
relatedness2_df <- relatedness2_df[
  with(relatedness2_df, order(relatedness2_df$unique_id)),
  ]
relatedness2_df <- subset(relatedness2_df, select = -c(INDV1, INDV2))
toDelete <- seq(1, nrow(relatedness2_df), 2)
relatedness2_df <- relatedness2_df[ toDelete ,]
```
####IBDMLE
```
ibd_mle_df <- ibd_mle_df %>%
  subset(INDV1 != INDV2) %>%
  select(-c(FileName)) 
ibd_mle_df <- ibd_mle_df[
  with(ibd_mle_df, order(ibd_mle_df$unique_id)),
  ]
ibd_mle_df <- subset(ibd_mle_df, select = -c(INDV1, INDV2))
toDelete <- seq(1, nrow(ibd_mle_df), 2)
ibd_mle_df <- ibd_mle_df[ toDelete ,]                  
```                         

# KING and IBS generation
```
gds_amur = snpgdsOpen("amurs-highcov-nofilter-biallelic-AN-MM-pcc.gds")
gds_bengal = snpgdsOpen("bengals-highcov-nofilter-biallelic-AN-MM-pcc.gds")
gds_malayan = snpgdsOpen("malayan-highcov-nofilter-biallelic-AN-MM-pcc.gds")
gds_sumatran = snpgdsOpen("sumatran-highcov-nofilter-biallelic-AN-MM-pcc.gds")
gds_indochinese = snpgdsOpen("indochinese-highcov-nofilter-biallelic-AN-MM-pcc.gds")
gds_generic = snpgdsOpen("generic-highcov-nofilter-biallelic-AN-MM-pcc.gds")
```
# Amur IBS/KING 
```
amur_sampIds = read.gdsn(index.gdsn(gds_amur, "sample.id")) #grab sample ids 
amur_KING <- snpgdsIBDKING(gds_amur, sample.id = amur_sampIds, autosome.only = F, num.thread = 4)
amur_KING_Matrix = amur_KING$IBS0
colnames(amur_KING_Matrix) <- rownames(amur_KING_Matrix) <- amur_KING$sample.id
amur_KING.df <- as.data.frame(as.table(amur_KING_Matrix))
amur_KING.df <- amur_KING.df %>%
  rowwise() %>%
  mutate(unique_id = paste(sort(c(as.character(Var1), as.character(Var2))), collapse = '-')) %>%
  ungroup() %>%
  subset(Var1 != Var2) %>%
  rename(KING = Freq)
amur_KING.df <- amur_KING.df[
  with(amur_KING.df, order(amur_KING.df$unique_id)),
  ]
amur_KING.df <- subset(amur_KING.df, select = -c(Var1, Var2))
toDelete <- seq(1, nrow(amur_KING.df), 2)
amur_KING.df <- amur_KING.df[ toDelete ,]   

amur_IBS = snpgdsIBS(gds_amur, sample.id = amur_sampIds, autosome.only = F, num.thread=4)
amur_IBS_Matrix = amur_IBS$ibs
colnames(amur_IBS_Matrix) <- rownames(amur_IBS_Matrix) <- amur_IBS$sample.id
amur_IBS.df <- as.data.frame(as.table(amur_IBS_Matrix ))
amur_IBS.df <- amur_IBS.df %>%
  rowwise() %>%
  mutate(unique_id = paste(sort(c(as.character(Var1), as.character(Var2))), collapse = '-')) %>%
  ungroup() %>%
  subset(Var1 != Var2) %>%
  rename(IBS = Freq)
amur_IBS.df <- amur_IBS.df[
  with(amur_IBS.df, order(amur_IBS.df$unique_id)),
  ]
amur_IBS.df <- subset(amur_IBS.df, select = -c(Var1, Var2))
toDelete <- seq(1, nrow(amur_IBS.df), 2)
amur_IBS.df <- amur_IBS.df[ toDelete ,]  
```
# Bengal IBS/KING 
```
bengal_sampIds = read.gdsn(index.gdsn(gds_bengal, "sample.id")) #grab sample ids 
bengal_KING <- snpgdsIBDKING(gds_bengal, sample.id = bengal_sampIds, autosome.only = F, num.thread = 4)
bengal_KING_Matrix = bengal_KING$IBS0
colnames(bengal_KING_Matrix) <- rownames(bengal_KING_Matrix) <- bengal_KING$sample.id
bengal_KING.df <- as.data.frame(as.table(bengal_KING_Matrix))
bengal_KING.df <- bengal_KING.df %>%
  rowwise() %>%
  mutate(unique_id = paste(sort(c(as.character(Var1), as.character(Var2))), collapse = '-')) %>%
  ungroup() %>%
  subset(Var1 != Var2) %>%
  rename(KING = Freq)
bengal_KING.df <- bengal_KING.df[
  with(bengal_KING.df, order(bengal_KING.df$unique_id)),
  ]
bengal_KING.df <- subset(bengal_KING.df, select = -c(Var1, Var2))
toDelete <- seq(1, nrow(bengal_KING.df), 2)
bengal_KING.df <- bengal_KING.df[ toDelete ,]   

bengal_IBS = snpgdsIBS(gds_bengal, sample.id = bengal_sampIds, autosome.only = F, num.thread=4)
bengal_IBS_Matrix = bengal_IBS$ibs
colnames(bengal_IBS_Matrix) <- rownames(bengal_IBS_Matrix) <- bengal_IBS$sample.id
bengal_IBS.df <- as.data.frame(as.table(bengal_IBS_Matrix ))
bengal_IBS.df <- bengal_IBS.df %>%
  rowwise() %>%
  mutate(unique_id = paste(sort(c(as.character(Var1), as.character(Var2))), collapse = '-')) %>%
  ungroup() %>%
  subset(Var1 != Var2) %>%
  rename(IBS = Freq)
bengal_IBS.df <- bengal_IBS.df[
  with(bengal_IBS.df, order(bengal_IBS.df$unique_id)),
  ]
bengal_IBS.df <- subset(bengal_IBS.df, select = -c(Var1, Var2))
toDelete <- seq(1, nrow(bengal_IBS.df), 2)
bengal_IBS.df <- bengal_IBS.df[ toDelete ,]  
```

# Malayan IBS/KING 
```
malayan_sampIds = read.gdsn(index.gdsn(gds_malayan, "sample.id")) #grab sample ids 
malayan_KING <- snpgdsIBDKING(gds_malayan, sample.id = malayan_sampIds, autosome.only = F, num.thread = 4)
malayan_KING_Matrix = malayan_KING$IBS0
colnames(malayan_KING_Matrix) <- rownames(malayan_KING_Matrix) <- malayan_KING$sample.id
malayan_KING.df <- as.data.frame(as.table(malayan_KING_Matrix))
malayan_KING.df <- malayan_KING.df %>%
  rowwise() %>%
  mutate(unique_id = paste(sort(c(as.character(Var1), as.character(Var2))), collapse = '-')) %>%
  ungroup() %>%
  subset(Var1 != Var2) %>%
  rename(KING = Freq)
malayan_KING.df <- malayan_KING.df[
  with(malayan_KING.df, order(malayan_KING.df$unique_id)),
  ]
malayan_KING.df <- subset(malayan_KING.df, select = -c(Var1, Var2))
toDelete <- seq(1, nrow(malayan_KING.df), 2)
malayan_KING.df <- malayan_KING.df[ toDelete ,]   

malayan_IBS = snpgdsIBS(gds_malayan, sample.id = malayan_sampIds, autosome.only = F, num.thread=4)
malayan_IBS_Matrix = malayan_IBS$ibs
colnames(malayan_IBS_Matrix) <- rownames(malayan_IBS_Matrix) <- malayan_IBS$sample.id
malayan_IBS.df <- as.data.frame(as.table(malayan_IBS_Matrix ))
malayan_IBS.df <- malayan_IBS.df %>%
  rowwise() %>%
  mutate(unique_id = paste(sort(c(as.character(Var1), as.character(Var2))), collapse = '-')) %>%
  ungroup() %>%
  subset(Var1 != Var2) %>%
  rename(IBS = Freq)
malayan_IBS.df <- malayan_IBS.df[
  with(malayan_IBS.df, order(malayan_IBS.df$unique_id)),
  ]
malayan_IBS.df <- subset(malayan_IBS.df, select = -c(Var1, Var2))
toDelete <- seq(1, nrow(malayan_IBS.df), 2)
malayan_IBS.df <- malayan_IBS.df[ toDelete ,]  
```

# Sumatran IBS/KING 
```
sumatran_sampIds = read.gdsn(index.gdsn(gds_sumatran, "sample.id")) #grab sample ids 
sumatran_KING <- snpgdsIBDKING(gds_sumatran, sample.id = sumatran_sampIds, autosome.only = F, num.thread = 4)
sumatran_KING_Matrix = sumatran_KING$IBS0
colnames(sumatran_KING_Matrix) <- rownames(sumatran_KING_Matrix) <- sumatran_KING$sample.id
sumatran_KING.df <- as.data.frame(as.table(sumatran_KING_Matrix))
sumatran_KING.df <- sumatran_KING.df %>%
  rowwise() %>%
  mutate(unique_id = paste(sort(c(as.character(Var1), as.character(Var2))), collapse = '-')) %>%
  ungroup() %>%
  subset(Var1 != Var2) %>%
  rename(KING = Freq)
sumatran_KING.df <- sumatran_KING.df[
  with(sumatran_KING.df, order(sumatran_KING.df$unique_id)),
  ]
sumatran_KING.df <- subset(sumatran_KING.df, select = -c(Var1, Var2))
toDelete <- seq(1, nrow(sumatran_KING.df), 2)
sumatran_KING.df <- sumatran_KING.df[ toDelete ,]   

sumatran_IBS = snpgdsIBS(gds_sumatran, sample.id = sumatran_sampIds, autosome.only = F, num.thread=4)
sumatran_IBS_Matrix = sumatran_IBS$ibs
colnames(sumatran_IBS_Matrix) <- rownames(sumatran_IBS_Matrix) <- sumatran_IBS$sample.id
sumatran_IBS.df <- as.data.frame(as.table(sumatran_IBS_Matrix ))
sumatran_IBS.df <- sumatran_IBS.df %>%
  rowwise() %>%
  mutate(unique_id = paste(sort(c(as.character(Var1), as.character(Var2))), collapse = '-')) %>%
  ungroup() %>%
  subset(Var1 != Var2) %>%
  rename(IBS = Freq)
sumatran_IBS.df <- sumatran_IBS.df[
  with(sumatran_IBS.df, order(sumatran_IBS.df$unique_id)),
  ]
sumatran_IBS.df <- subset(sumatran_IBS.df, select = -c(Var1, Var2))
toDelete <- seq(1, nrow(sumatran_IBS.df), 2)
sumatran_IBS.df <- sumatran_IBS.df[ toDelete ,]  
```
# Indochinese IBS/KING 
```
indochinese_sampIds = read.gdsn(index.gdsn(gds_indochinese, "sample.id")) #grab sample ids 
indochinese_KING <- snpgdsIBDKING(gds_indochinese, sample.id = indochinese_sampIds, autosome.only = F, num.thread = 4)
indochinese_KING_Matrix = indochinese_KING$IBS0
colnames(indochinese_KING_Matrix) <- rownames(indochinese_KING_Matrix) <- indochinese_KING$sample.id
indochinese_KING.df <- as.data.frame(as.table(indochinese_KING_Matrix))
indochinese_KING.df <- indochinese_KING.df %>%
  rowwise() %>%
  mutate(unique_id = paste(sort(c(as.character(Var1), as.character(Var2))), collapse = '-')) %>%
  ungroup() %>%
  subset(Var1 != Var2) %>%
  rename(KING = Freq)
indochinese_KING.df <- indochinese_KING.df[
  with(indochinese_KING.df, order(indochinese_KING.df$unique_id)),
  ]
indochinese_KING.df <- subset(indochinese_KING.df, select = -c(Var1, Var2))
toDelete <- seq(1, nrow(indochinese_KING.df), 2)
indochinese_KING.df <- indochinese_KING.df[ toDelete ,]   

indochinese_IBS = snpgdsIBS(gds_indochinese, sample.id = indochinese_sampIds, autosome.only = F, num.thread=4)
indochinese_IBS_Matrix = indochinese_IBS$ibs
colnames(indochinese_IBS_Matrix) <- rownames(indochinese_IBS_Matrix) <- indochinese_IBS$sample.id
indochinese_IBS.df <- as.data.frame(as.table(indochinese_IBS_Matrix ))
indochinese_IBS.df <- indochinese_IBS.df %>%
  rowwise() %>%
  mutate(unique_id = paste(sort(c(as.character(Var1), as.character(Var2))), collapse = '-')) %>%
  ungroup() %>%
  subset(Var1 != Var2) %>%
  rename(IBS = Freq)
indochinese_IBS.df <- indochinese_IBS.df[
  with(indochinese_IBS.df, order(indochinese_IBS.df$unique_id)),
  ]
indochinese_IBS.df <- subset(indochinese_IBS.df, select = -c(Var1, Var2))
toDelete <- seq(1, nrow(indochinese_IBS.df), 2)
indochinese_IBS.df <- indochinese_IBS.df[ toDelete ,]  
```
# Generic IBS/KING 
```
generic_sampIds = read.gdsn(index.gdsn(gds_generic, "sample.id")) #grab sample ids 
generic_KING <- snpgdsIBDKING(gds_generic, sample.id = generic_sampIds, autosome.only = F, num.thread = 4)
generic_KING_Matrix = generic_KING$IBS0
colnames(generic_KING_Matrix) <- rownames(generic_KING_Matrix) <- generic_KING$sample.id
generic_KING.df <- as.data.frame(as.table(generic_KING_Matrix))
generic_KING.df <- generic_KING.df %>%
  rowwise() %>%
  mutate(unique_id = paste(sort(c(as.character(Var1), as.character(Var2))), collapse = '-')) %>%
  ungroup() %>%
  subset(Var1 != Var2) %>%
  rename(KING = Freq)
generic_KING.df <- generic_KING.df[
  with(generic_KING.df, order(generic_KING.df$unique_id)),
  ]
generic_KING.df <- subset(generic_KING.df, select = -c(Var1, Var2))
toDelete <- seq(1, nrow(generic_KING.df), 2)
generic_KING.df <- generic_KING.df[ toDelete ,]   

generic_IBS = snpgdsIBS(gds_generic, sample.id = generic_sampIds, autosome.only = F, num.thread=4)
generic_IBS_Matrix = generic_IBS$ibs
colnames(generic_IBS_Matrix) <- rownames(generic_IBS_Matrix) <- generic_IBS$sample.id
generic_IBS.df <- as.data.frame(as.table(generic_IBS_Matrix ))
generic_IBS.df <- generic_IBS.df %>%
  rowwise() %>%
  mutate(unique_id = paste(sort(c(as.character(Var1), as.character(Var2))), collapse = '-')) %>%
  ungroup() %>%
  subset(Var1 != Var2) %>%
  rename(IBS = Freq)
generic_IBS.df <- generic_IBS.df[
  with(generic_IBS.df, order(generic_IBS.df$unique_id)),
  ]
generic_IBS.df <- subset(generic_IBS.df, select = -c(Var1, Var2))
toDelete <- seq(1, nrow(generic_IBS.df), 2)
generic_IBS.df <- generic_IBS.df[ toDelete ,]  

IBS_all <- do.call("rbind", list(amur_IBS.df, malayan_IBS.df, sumatran_IBS.df, indochinese_IBS.df, bengal_IBS.df, generic_IBS.df))
KING_all <- do.call("rbind", list(amur_KING.df, malayan_KING.df, sumatran_KING.df, indochinese_KING.df, bengal_KING.df, generic_KING.df))
```
# Merging 
```
vcftools <- merge(relatedness2_df,relatedness_df, by="unique_id")
IBS_KING <- merge(IBS_all,KING_all, by="unique_id")
vcftools_mle <- merge(vcftools, ibd_mle_df, by="unique_id")  
vcftools_mle_IBS_KING <- merge(vcftools_mle, IBS_KING, by="unique_id")
vcftools_mle_IBS_KING_pedigree <- full_join(vcftools_mle_IBS_KING, amur_mal_sum_pedigree, by="unique_id")                    
```
#transform data to 0's
```
vcftools_mle_IBS_KING_pedigree[, c(2)][vcftools_mle_IBS_KING_pedigree[, c(2)] < 0] <- 0
vcftools_mle_IBS_KING_pedigree[, c(4)][vcftools_mle_IBS_KING_pedigree[, c(4)] < 0] <- 0
```
#sort for easy plotting
```
vcftools_mle_IBS_KING_pedigree <- vcftools_mle_IBS_KING_pedigree[, c(1, 3, 2, 4, 5, 6, 7, 8)]
```
#First plot correlations 
```
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Sumatran" = "cornflowerblue", "Indochinese" = "gold4", "Generic"="gray25")#palette
```

# Final correlation plot 
#rename columns for correlations plot
```
vcftools_mle_IBS_KING_pedigree <- vcftools_mle_IBS_KING_pedigree %>%
  rename(Vcftools2 = RELATEDNESS_PHI, Vcftools1 = RELATEDNESS_AJK, 'IBD MLE SNPRelate' = GDS_MLE_Kinship, 'Pedigree Kinship' = pedigree_kinship)

library(GGally)
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_manual(name = "Subspecies", values = cbPalette) + scale_fill_manual(name = "Subspecies", values = cbPalette)
unlockBinding("ggplot",parent.env(asNamespace("GGally")))
assign("ggplot",ggplot,parent.env(asNamespace("GGally")))

ggpairs(vcftools_mle_IBS_KING_pedigree, columns = 3:8, 
        ggplot2::aes(colour = Subspecies))
```
#now plot IBD MLE individuals
#first sort according to subspecies so it looks pretty
```
x <- c('Amur', 'Bengal', 'Indochinese', 'Malayan', 'Sumatran', 'Generic')
vcftools_mle_IBS_KING_pedigree$Subspecies <- factor(vcftools_mle_IBS_KING_pedigree$Subspecies, 
                                                levels=c('Amur', 'Bengal', 'Indochinese', 'Malayan', 'Sumatran', 'Generic'))


vcftools_mle_IBS_KING_pedigree %>%
  select(c(unique_id,Subspecies,GDS_MLE_Kinship)) %>%
  arrange(Subspecies) %>%
  mutate(unique_id = forcats::fct_inorder(factor(unique_id))) %>%
  ggplot(aes(x=unique_id, y=GDS_MLE_Kinship, col = Subspecies)) +
  geom_point() +
  scale_color_manual(values = cbPalette) +
  geom_hline(yintercept = 0.5, linetype = "solid", color = "goldenrod", size = 1) +
  geom_hline(yintercept = 0.354, linetype = "dashed", color = "goldenrod", size = 1) +
  geom_hline(yintercept = 0.25, linetype = "solid", color = "lightskyblue", size = 1) +
  geom_hline(yintercept = 0.177, linetype = "dashed", color = "lightskyblue", size = 1) +
  geom_hline(yintercept = 0.125, linetype = "solid", color = "mediumvioletred", size = 1) +
  geom_hline(yintercept = 0.0884, linetype = "dashed", color = "mediumvioletred", size = 1) +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab('IBD MLE Kinship (SNPRelate)') + 
  xlab('Pairs of Individuals') 
```
