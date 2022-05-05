#Load Libraries
library(tidyverse)
library(data.table)
library(ggpubr)
options(scipen = 999)

###Function to make a folded SFS###
#THIS IS BERNARD's ORIGINAL >> need to add 0 to front and end to stand in for monomorphic (don't need to do this if it came from dadi)
fold <- function(SFSCountCol, n, norm=TRUE){
  if (length(SFSCountCol) < (n+1)){
    data = c(SFSCountCol, rep(0,(n+1)-length(SFSCountCol)))
  }
  data = SFSCountCol[2:n] # adds together sfs and backwards sfs
  data_fold = data + rev(data) # takes first half of that added together sfs (but not middle entry)
  data_fold = data_fold[1:(n/2-1)]# adds middle entry that didn't have anything added to the end
  data_fold = c(data_fold,data[(n/2)])# truncates and sums up remaining fields if desired (not needed here)
  #data_trunc = c(data_fold[1:(trunc-1)],sum(data_fold[trunc:length(data_fold)]))
  #if (norm){
  #  data_trunc = data_trunc/sum(data_trunc)
  #}
  #return(data_trunc)
  return(tibble(data_fold))
}

####Start with N10
setwd("~/Documents/Tigers/sfs/unFolded/OnlySIFTAnnots/")
fnames = list.files(pattern = "\\N10.txt$")
fnames_N6 = list.files(pattern = "\\N6.txt$")

unfoldedDF = rbindlist(sapply(fnames, read.delim, simplify = FALSE), use.names = TRUE, idcol = "Subspecies") %>%
  mutate(Subspecies = gsub("_.*$","",Subspecies),
         Annot = ifelse(is.na(FreqBin_SY), "Putatively deleterious", "Putatively neutral"),
         FreqBin = ifelse(is.na(FreqBin_SY), FreqBin_Del,FreqBin_SY),
         FreqBin_SY = NULL,
         FreqBin_Del = NULL) 

unfoldedDF_N6 = rbindlist(sapply(fnames_N6, read.delim, simplify = FALSE), use.names = TRUE, idcol = "Subspecies") %>%
  mutate(Subspecies = gsub("_.*$","",Subspecies),
         Annot = ifelse(is.na(FreqBin_SY), "Putatively deleterious", "Putatively neutral"),
         FreqBin = ifelse(is.na(FreqBin_SY), FreqBin_Del,FreqBin_SY),
         FreqBin_SY = NULL,
         FreqBin_Del = NULL) 

#make a data frame with all spots filled
Subspecies = rep(c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'Sumatran'), each = 21)
full = Subspecies %>%
  as.data.frame() %>%
  mutate(FreqBin = rep(seq(from = 0, to = 20, by = 1), length.out = n())) 
colnames(full) = c("Subspecies", "FreqBin")

Subspecies = rep(c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'Sumatran'), each = 13)
full_N6 = Subspecies %>%
  as.data.frame() %>%
  mutate(FreqBin = rep(seq(from = 0, to = 12, by = 1), length.out = n())) 
colnames(full_N6) = c("Subspecies", "FreqBin")

#Reformatting
PlottingUnfolded = unfoldedDF %>%
  filter(FreqBin > 0 & FreqBin < 20) %>% #remove fixed stuff
  group_by(Subspecies, FreqBin, Annot) %>%
  summarise(WholeGenomeCounts = sum(sum)) %>%
  ungroup() %>%
  mutate(Subspecies = factor(Subspecies, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'Sumatran')))

TotalCounts = PlottingUnfolded %>%
  group_by(Subspecies, Annot) %>%
  summarise(TotalSites = sum(WholeGenomeCounts))

PlottingUnfolded = PlottingUnfolded %>%
  inner_join(TotalCounts) %>%
  mutate(Proportional = WholeGenomeCounts/TotalSites)

###Fold it
split_N10 = unfoldedDF %>%
  group_split(Annot)

PlottingFolded_del = split_N10[[1]] %>%
  group_by(Subspecies, FreqBin) %>%
  summarise(WholeGenomeCounts = sum(sum)) %>%
  ungroup() %>%
  right_join(full) %>%
  filter(Subspecies != "Indochinese") %>% 
  replace(is.na(.), as.integer(0)) %>%
  arrange(Subspecies, FreqBin) %>% #fxn needs data to be sorted from 0:2n for each subspecies
  group_by(Subspecies) %>%
  group_modify(~ fold(.x$WholeGenomeCounts, n = 20)) %>% #this n = num chroms
  ungroup() %>%
  mutate(bin = rep(seq(1:10), length.out = n()),
         Subspecies = factor(Subspecies, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'Sumatran')),
         Annot = "Putatively deleterious")

PlottingFolded = split_N10[[2]] %>%
  group_by(Subspecies, FreqBin) %>%
  summarise(WholeGenomeCounts = sum(sum)) %>%
  ungroup() %>%
  right_join(full) %>%
  filter(Subspecies != "Indochinese") %>% 
  replace(is.na(.), 0) %>%
  arrange(Subspecies, FreqBin) %>% #fxn needs data to be sorted from 0:2n for each subspecies
  group_by(Subspecies) %>%
  group_modify(~ fold(.x$WholeGenomeCounts, n = 20)) %>% #this n = num chroms
  ungroup() %>%
  mutate(bin = rep(seq(1:10), length.out = n()),
         Subspecies = factor(Subspecies, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'Sumatran')),
         Annot = "Putatively neutral") %>%
  rbind.data.frame(PlottingFolded_del) %>%
  inner_join(TotalCounts) %>%
  mutate(Proportional = data_fold/TotalSites)

###N6
PlottingUnfolded_N6 = unfoldedDF_N6 %>%
  filter(FreqBin > 0 & FreqBin < 12) %>% #remove fixed stuff
  group_by(Subspecies, FreqBin, Annot) %>%
  summarise(WholeGenomeCounts = sum(sum)) %>%
  ungroup() %>%
  mutate(Subspecies = factor(Subspecies, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'Sumatran')))

TotalCounts_N6 = PlottingUnfolded_N6 %>%
  group_by(Subspecies, Annot) %>%
  summarise(TotalSites = sum(WholeGenomeCounts))

PlottingUnfolded_N6 = PlottingUnfolded_N6 %>%
  inner_join(TotalCounts_N6) %>%
  mutate(Proportional = WholeGenomeCounts/TotalSites)

###Fold it
split_N6 = unfoldedDF_N6 %>%
  group_split(Annot)

PlottingFolded_del_N6 = split_N6[[1]] %>%
  group_by(Subspecies, FreqBin) %>%
  summarise(WholeGenomeCounts = sum(sum)) %>%
  ungroup() %>%
  right_join(full_N6) %>%
  replace(is.na(.), as.integer(0)) %>%
  arrange(Subspecies, FreqBin) %>% #fxn needs data to be sorted from 0:2n for each subspecies
  group_by(Subspecies) %>%
  group_modify(~ fold(.x$WholeGenomeCounts, n = 12)) %>% #this n = num chroms
  ungroup() %>%
  mutate(bin = rep(seq(1:6), length.out = n()),
         Subspecies = factor(Subspecies, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'Sumatran')),
         Annot = "Putatively deleterious")

PlottingFolded_N6 = split_N6[[2]] %>%
  group_by(Subspecies, FreqBin) %>%
  summarise(WholeGenomeCounts = sum(sum)) %>%
  ungroup() %>%
  right_join(full_N6) %>%
  replace(is.na(.), 0) %>%
  arrange(Subspecies, FreqBin) %>% #fxn needs data to be sorted from 0:2n for each subspecies
  group_by(Subspecies) %>%
  group_modify(~ fold(.x$WholeGenomeCounts, n = 12)) %>% #this n = num chroms
  ungroup() %>%
  mutate(bin = rep(seq(1:6), length.out = n()),
         Subspecies = factor(Subspecies, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'Sumatran')),
         Annot = "Putatively neutral") %>%
  rbind.data.frame(PlottingFolded_del_N6) %>%
  inner_join(TotalCounts_N6) %>%
  mutate(Proportional = data_fold/TotalSites)



####Plotting
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Sumatran" = "cornflowerblue", "Indochinese" = "gold4", "Generic"="gray25")#palette

UnfoldedSFS_N10 = ggplot(PlottingUnfolded, aes(y=Proportional, x=FreqBin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_x_continuous(breaks=1:19) + 
  facet_grid(Annot~.) +
  ylim(0,0.3)+
  scale_fill_manual(values= cbPalette) + 
  labs(x = "SNP Frequency", 
       y= "Proportion of SNPs", 
       title = "N=10") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20),
        strip.text = element_text(size  = 14),
        plot.title=element_text(size  = 24, hjust = 0.5), 
        axis.title=element_text(size  = 24),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "none")

FoldedSFS_N10 = ggplot(PlottingFolded, aes(y=Proportional, x=bin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  facet_grid(Annot~., scales = "free_y") +
  scale_x_continuous(breaks=1:10) + 
  ylim(0,0.3)+
  scale_fill_manual(values= cbPalette) + 
  labs(x = "SNP Frequency", 
       y= "Counts of SNPs", 
       title = "Folded SFS (based on SIFT annotations)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20),
        strip.text = element_text(size  = 14),
        plot.title=element_text(size  = 24, hjust = 0.5), 
        axis.title=element_text(size  = 24),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "none") 

UnfoldedSFS_N6 = ggplot(PlottingUnfolded_N6, aes(y=Proportional, x=FreqBin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_x_continuous(breaks=1:11) + 
  facet_grid(Annot~.) +
  ylim(0,0.3)+
  scale_fill_manual(values= cbPalette) + 
  labs(x = "SNP Frequency", 
       y= "Proportion of SNPs", 
       title = "N=6") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20),
        strip.text = element_text(size  = 14),
        plot.title=element_text(size  = 24, face = "bold", hjust = 0.5), 
        axis.title=element_text(size  = 24),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = c(0.8, 0.8)) 

FoldedSFS_N6 = ggplot(PlottingFolded_N6, aes(y=Proportional, x=bin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_x_continuous(breaks=1:6) + 
  facet_grid(Annot~.) +
  ylim(0,0.3) +
  scale_fill_manual(values= cbPalette) + 
  labs(x = "SNP Frequency", 
       y= "Proportion of SNPs", 
       title = "All Populations Whole Genome SFS") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20),
        strip.text = element_text(size  = 14),
        plot.title=element_text(size  = 24, hjust = 0.5), 
        axis.title=element_text(size  = 24),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "none") 


allSFS = ggarrange(UnfoldedSFS_N10, UnfoldedSFS_N6, FoldedSFS_N10, FoldedSFS_N6, 
                   labels = c("A","B","C","D"), font.label = list(size = 18), 
                   nrow = 2, ncol = 2)

allSFS_axes = annotate_figure(allSFS, 
                                  left = text_grob("Proportion of variable sites", 
                                                   color = "black", 
                                                   face = "bold", 
                                                   size = 24, 
                                                   rot = 90), 
                                  bottom = text_grob("Site frequency", 
                                                     color = "black", 
                                                     face = "bold", 
                                                     size = 24))

#jpeg(file="~/Documents/Tigers/sfs/SuppSFS.jpeg", width = 2200, height = 1600, res = 100)
print(allSFS_axes)
#dev.off()
