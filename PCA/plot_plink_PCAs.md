---
Plot PCAs in Plink
---

Load libraries

```
library(ggpubr)
library(tidyverse)
```
Load in data files from plink and subspecies assignments
```
eigenvec_highcov_nocorr <- read.table('highcov-nofilter-biallelic-AN-MM-pcc.eigenvec',head=F)
eigenval_highcov_nocorr <- read.table('highcov-nofilter-biallelic-AN-MM-pcc.eigenval',head=F)
metadata <- read.table('~/Documents/Documents - Ellie’s MacBook Pro (2)/captives-new-2021/Metadata-files/plotting_metadata.csv', 
                       header=TRUE, sep = ',')
```
Set colors for plot
```
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", 
              "Sumatran" = "cornflowerblue", "Indochinese" = "gold4", 
              "South China" = "plum", "Generic"="gray25")#palette
              
```
Calculate percent variance from eigen values. This should be the eigenvalue divided by the number of loadings (individuals) *100 for percentage
```
eigenval_highcov_nocorr$V2 <- (eigenval_highcov_nocorr$V1/171)*100
```

Clean up data (rename columns so they match and merge with metadata for plotting)
```
eigenvec_table_rename <- eigenvec_highcov_nocorr %>% 
  rename(Individual=V2)
eigenvec_with_sub <- merge(eigenvec_table_rename, metadata,by="Individual")

```
Plot and set output
```
PC1 <- ggplot(eigenvec_with_sub, aes(x=V3, y=V4, color=Subspecies_GroupID)) +
  scale_color_manual(name="Subspecies", values=cbPalette) +
  geom_point() +
  xlab("PC1 (7.13%)") + ylab("PC2 (6.22%)") +
  theme(legend.position = "right", legend.title = element_blank()) +
  theme_bw()

PC2 <- ggplot(eigenvec_with_sub, aes(x=V3, y=V5, color=Subspecies_GroupID)) +
  scale_color_manual(name="Subspecies", values=cbPalette) +
  geom_point() +
  xlab("PC1 (7.13%)") + ylab("PC3 (5.35%)") +
  theme_bw() +
  theme(legend.position = "none") 

PC3 <- ggplot(eigenvec_with_sub, aes(x=V4, y=V5, color=Subspecies_GroupID)) +
  scale_color_manual(name="Subspecies", values=cbPalette) +
  geom_point() +
  xlab("PC2 (6.22%)") + ylab("PC3 (5.35%)") +
  theme_bw() +
  theme(legend.position = "none")


figure <- ggarrange(PC1, ggarrange(PC2, PC3, ncol=2, labels = c("B" , "C")),
                    nrow = 2, labels = "A")

jpeg(file="~/Documents/Documents - Ellie’s MacBook Pro (2)/captives-new-2021/pca/pca-highcov-nodups.jpeg", width = 900, height = 700, res = 100)
```
