#Load Libraries
library(tidyverse)
library(data.table)
library(ggpubr)
options(scipen = 999)
setwd("~/Documents/Tigers/sfs/unFolded/")

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
fnames = list.files(pattern = "\\N10.txt$")
fnames_N6 = list.files(pattern = "\\N6.txt$")

unfoldedDF = rbindlist(sapply(fnames, read.delim, simplify = FALSE), use.names = TRUE, idcol = "Subspecies") %>%
  mutate(Subspecies = gsub("_.*$","",Subspecies))

unfoldedDF_N6 = rbindlist(sapply(fnames_N6, read.delim, simplify = FALSE), use.names = TRUE, idcol = "Subspecies") %>%
  mutate(Subspecies = gsub("_.*$","",Subspecies))

###N6
#run when group is missing bins
# Indochinese = unfoldedDF_N6 %>%
#   filter(Subspecies == "Indochinese") 
# 
# Indochinese_addBins = unfoldedDF_N6 %>%
#   filter(Subspecies == "Amur") %>%
#   mutate(sum = 0,
#          Subspecies = "Indochinese") %>%
#   left_join(Indochinese, by = c("Subspecies", "FreqBin","chromosome")) %>%
#   mutate(sum = ifelse(is.na(sum.y), sum.x, sum.y)) %>%
#   select("Subspecies", "FreqBin", "sum", "chromosome")
  

PlottingFolded_N6 = unfoldedDF_N6 %>%
  #filter(Subspecies != "Indochinese") %>%
  #rbind.data.frame(Indochinese_addBins) %>%
  group_by(Subspecies, FreqBin) %>% 
  summarise(WholeGenomeCounts = sum(sum)) %>%
  ungroup() %>%
  group_by(Subspecies) %>%
  group_modify(~ fold(.x$WholeGenomeCounts, n = 12)) %>% #this n = num chroms
  ungroup() %>%
  mutate(bin = rep(seq(1:6), length.out = n()),
         Subspecies = factor(Subspecies, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'Sumatran')))

PlottingUnfolded_N6 = unfoldedDF_N6 %>%
  #filter(Subspecies != "Indochinese") %>%
  #rbind.data.frame(Indochinese_addBins) %>%
  filter(FreqBin > 0 & FreqBin < 12) %>% #remove fixed stuff
  group_by(Subspecies, FreqBin) %>%
  summarise(WholeGenomeCounts = sum(sum)) %>%
  ungroup() %>%
  mutate(Subspecies = factor(Subspecies, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'Sumatran')))

TotalCounts_N6 = PlottingUnfolded_N6 %>%
  group_by(Subspecies) %>%
  summarise(TotalSites = sum(WholeGenomeCounts))

PlottingUnfolded_N6$Proportional = PlottingUnfolded_N6$WholeGenomeCounts/(TotalCounts_N6$TotalSites[match(PlottingUnfolded_N6$Subspecies, TotalCounts_N6$Subspecies)])

PlottingFolded_N6$Proportional = PlottingFolded_N6$data_fold/(TotalCounts_N6$TotalSites[match(PlottingFolded_N6$Subspecies, TotalCounts_N6$Subspecies)])

####Plotting
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Sumatran" = "cornflowerblue", "Indochinese" = "gold4", "Generic"="gray25")#palette

UnfoldedSFS_N10 = ggplot(PlottingUnfolded, aes(y=Proportional, x=FreqBin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_x_continuous(breaks=1:19) + 
  ylim(0,0.5)+
  scale_fill_manual(values= cbPalette) + 
  labs(x = "SNP Frequency", 
       y= "Proportion of SNPs", 
       title = "N=10") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        plot.title=element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title=element_blank(),
        legend.position = "none") 

FoldedSFS_N10 = ggplot(PlottingFolded, aes(y=Proportional, x=bin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_x_continuous(breaks=1:10) + 
  ylim(0,0.5)+
  scale_fill_manual(values= cbPalette) + 
  labs(x = "SNP Frequency", 
       y= "Proportion of SNPs", 
       title = "All Populations Whole Genome SFS") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        plot.title=element_blank(), 
        axis.title=element_blank(),
        legend.position = "none") 

UnfoldedSFS_N6 = ggplot(PlottingUnfolded_N6, aes(y=Proportional, x=FreqBin, fill=Subspecies)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_x_continuous(breaks=1:11) + 
  ylim(0,0.5)+
  scale_fill_manual(values= cbPalette) + 
  labs(x = "SNP Frequency", 
       y= "Proportion of SNPs", 
       title = "N=6") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        plot.title=element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title=element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.position = c(0.8, 0.8)) 

FoldedSFS_N6 = ggplot(PlottingFolded_N6, aes(y=Proportional, x=bin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_x_continuous(breaks=1:6) + 
  ylim(0,0.5)+
  scale_fill_manual(values= cbPalette) + 
  labs(x = "SNP Frequency", 
       y= "Proportion of SNPs", 
       title = "All Populations Whole Genome SFS") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        plot.title=element_blank(), 
        axis.title=element_blank(),
        legend.position = "none") 


allSFS = ggarrange(UnfoldedSFS_N10, UnfoldedSFS_N6, FoldedSFS_N10, FoldedSFS_N6, 
                   #labels = c("A","B","C","D"), font.label = list(size = 18), 
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
