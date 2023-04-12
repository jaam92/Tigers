---
Plot PCAs in Plink
---

Load libraries

```
library(ggpubr)
library(tidyverse)
library(patchwork)
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

# geography and new ancestry pca ------------------------------------------
generic_geo_vec <- read.table('~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Ellie’s MacBook Pro (2)/Captives-2022-Final/PCA/rescue_pca_highlowcov.eigenvec',head=F)
generic_geo_val <- read.table('~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Ellie’s MacBook Pro (2)/Captives-2022-Final/PCA/rescue_pca_highlowcov.eigenval',head=F)

generic_geo_val$V2 <- (generic_geo_val$V1/113)*100

generic_geo_vec <- generic_geo_vec %>% 
  rename(Individual=V2)
eigenvec_geo <- merge(generic_geo_vec, metadata,by="Individual")
library(randomcoloR)
n <- 23
palette <- distinctColorPalette(n)

a<- ggplot(eigenvec_geo, aes(x=V3, y=V4, color=Top_Ancestry)) +
  scale_color_manual(name="Top Ancestry", values=Ancestry_palette) +
  geom_point(size = 2) +
  xlab("PC1 (3.04%)") + ylab("PC2 (2.35%)") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.25), axis.text.x = element_text(angle = 60, hjust=1, size = 12), 
        axis.text.y = element_text(size = 12), axis.title.x=element_text(size = 12), axis.title.y = element_text(size= 12)) 


b <- ggplot(eigenvec_geo, aes(x=V4, y=V5, color=Top_Ancestry)) +
  scale_color_manual(name="Rescue_Location", values=Ancestry_palette) +
  geom_point(aes(shape=Top_Ancestry), size = 2) +
  xlab("PC2 (2.36%)") + ylab("PC3 (2.26%)") +
  theme(legend.position = "right", legend.title = "Top Ancestry") +
  theme_bw() +
  theme(legend.position = "none") 

n <- 6
palette2 <- distinctColorPalette(n)


c<- ggplot(eigenvec_geo, aes(x=V3, y=V4, color=Region)) +
  scale_color_brewer(palette = "Spectral", name="Region") +
  geom_point(size = 2) +
  xlab("PC1 (3.04%)") + ylab("PC2 (2.35%)") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.25), axis.text.x = element_text(angle = 60, hjust=1, size = 12), 
        axis.text.y = element_text(size = 12), axis.title.x=element_text(size = 12), axis.title.y = element_text(size= 12)) 

d <- ggplot(eigenvec_geo, aes(x=V4, y=V5, color=Region)) +
  scale_color_manual(name="Region", values=palette2) +
  geom_point(aes(shape=Top_Ancestry), size = 2) +
  xlab("PC2 (2.36%)") + ylab("PC3 (2.26%)") +
  theme(legend.position = "right", legend.title = "element_blank()") +
  theme_bw() +
  theme(legend.position = "none") 
a|c

write.csv(eigenvec_geo, '~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Ellie’s MacBook Pro (2)/Captives-2022-Final/PCA/eigenvec_geo.csv')
```
