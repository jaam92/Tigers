#Load Libraries
library(tidyverse)
library(data.table)
options(scipen = 999)
setwd("~/Documents/Tigers/sfs")

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

####unfolded sfs
fnames = list.files(pattern = "\\N10.txt$")

unfoldedDF = rbindlist(sapply(fnames, read.delim, simplify = FALSE), use.names = TRUE, idcol = "Subspecies") %>%
  mutate(Subspecies = gsub("_.*$","",Subspecies))

#Plot Data
PlottingFolded= unfoldedDF %>%
  group_by(Subspecies, FreqBin) %>%
  summarise(WholeGenomeCounts = sum(sum)) %>%
  ungroup() %>%
  group_by(Subspecies) %>%
  group_modify(~ fold(.x$WholeGenomeCounts, n = 20)) %>% #this n = num chroms
  ungroup() %>%
  mutate(bin = rep(seq(1:10), length.out = n()))

PlottingUnfolded = unfoldedDF %>%
  filter(FreqBin > 0 & FreqBin < 20) %>% #remove fixed stuff
  group_by(Subspecies, FreqBin) %>%
  summarise(WholeGenomeCounts = sum(sum)) %>%
  ungroup() 

TotalCounts = PlottingUnfolded %>%
  group_by(Subspecies) %>%
  summarise(TotalSites = sum(WholeGenomeCounts))

PlottingUnfolded$Proportional = PlottingUnfolded$WholeGenomeCounts/(TotalCounts$TotalSites[match(PlottingUnfolded$Subspecies, TotalCounts$Subspecies)])

PlottingFolded$Proportional = PlottingFolded$data_fold/(TotalCounts$TotalSites[match(PlottingFolded$Subspecies, TotalCounts$Subspecies)])

#Plot only the bins starting with singletons
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Sumatran" = "cornflowerblue", "Indochinese" = "gold4", "South China" = "plum", "Generic"="gray25")#palette

UnfoldedSFS_N10 = ggplot(PlottingUnfolded, aes(y=Proportional, x=FreqBin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_x_continuous(breaks=1:19) + 
  scale_fill_manual(values= cbPalette) + 
  labs(x = "SNP Frequency", 
       y= "Proportion of SNPs", 
       title = "All Populations Whole Genome SFS") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) 

FoldedSFS_N10 = ggplot(PlottingFolded, aes(y=Proportional, x=bin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_x_continuous(breaks=1:10) + 
  scale_fill_manual(values= cbPalette) + 
  labs(x = "SNP Frequency", 
       y= "Proportion of SNPs", 
       title = "All Populations Whole Genome SFS") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) 
