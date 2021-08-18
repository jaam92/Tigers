#Load Libraries
library(tidyverse)
library(data.table)
setwd("~/Documents/TigerProject/Rscripts/sfs/Unfolded/")

#Load files
fnames = list.files(pattern = "\\_N18.txt$")

df = rbindlist(sapply(fnames, read.delim, simplify = FALSE), use.names = TRUE, idcol = "Subspecies") %>%
  mutate(Subspecies = gsub("_.*$","",Subspecies)) %>%
  group_by(Subspecies, FreqBin) %>%
  summarise(WholeGenomeCounts = sum(sum)) %>%
  ungroup()

fixedSites = read_delim("allPutCatChroms_fixedSites_polarized.bed", delim = "\t") %>%
  mutate(numBP = endVCF-startVCF,
         DerHOM = ifelse(TigerAllele == DER, 1, 0))

uniqueSites = fixedSites[!duplicated(fixedSites[,6:8]),] %>%
  group_by(DerHOM) %>%
  summarise_at(c("numBP"), sum, na.rm = TRUE)

duplicatedSites = fixedSites[duplicated(fixedSites[,6:8]),]

#reformat data
df$bin = rep(seq(0:36), length(fnames)) #add bin number

TotalCountsFolded = df %>%
  group_by(Subspecies) %>%
  summarise(TotalSites = sum(data_fold))

df$Proportional = df$data_fold/(TotalCountsFolded$TotalSites[match(df$Subspecies, TotalCountsFolded$Subspecies)]) #add column with proportion of sites

#Now plot
cbPalette = c("Generic" = "gray25", "Amur" = "#D55E00",  "Bengal" = "steelblue", "Malayan" = "#009E73", "Sumatran" = "gold3", "Indochinese" = "chocolate1", "Sanctuary" = "#CC79A7", "F1Wild" = "#867BCF", "Wild"="goldenrod2", "F2Wild"="coral2", "Captive"="cornflowerblue", "Unknown"="darkmagenta", "Zoo"="darkseagreen3") #palette


propFoldedSFS_withGeneric = ggplot(df, aes(y=Proportional, x=bin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_fill_manual(values= cbPalette, name = "Sub-Species") + 
  scale_x_continuous(breaks = df$bin) + 
  labs(x= "SNP Frequency", y= "Proportion of SNPs") + 
  ggtitle("All Sub-Species (N=18)\nWhole Genome Folded SFS") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 10),
        plot.title=element_text(size=24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18))

countsFoldedSFS_withGeneric = ggplot(df, aes(y=data_fold, x=bin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_fill_manual(values= cbPalette, name = "Sub-Species") + 
  scale_x_continuous(breaks = df$bin) + 
  labs(x= "SNP Frequency", y= "Proportion of SNPs") + 
  ggtitle("All Sub-Species (N=18)\nWhole Genome Folded SFS") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 10),
        plot.title=element_text(size=24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18))

propFoldedSFS_withoutGeneric = ggplot(df %>% 
                                        filter(Subspecies != "Generic"),
                                      aes(y=Proportional, x=bin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_fill_manual(values= cbPalette, name = "Sub-Species") + 
  scale_x_continuous(breaks = df$bin) + 
  labs(x= "SNP Frequency", y= "Proportion of SNPs") + 
  ggtitle("All Sub-Species (N=18)\nWhole Genome Folded SFS") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 10),
        plot.title=element_text(size=24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18))
