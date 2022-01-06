#Load Libraries
library(tidyverse)
library(data.table)
options(scipen = 999)
setwd("~/TigerProject/sfs/Folded/")

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

#Load files
fnames = list.files(pattern = "\\N18_indivs$")

df = rbindlist(sapply(fnames, read.delim, simplify = FALSE), use.names = TRUE, idcol = "Subspecies") %>%
  mutate(Subspecies = gsub("_.*$","",Subspecies)) %>%
  group_by(Subspecies) %>%
  group_modify(~ fold(.x$count, n = 18))
  

df$bin = rep(seq(1:9), length(fnames)) #add bin number

TotalCountsFolded = df %>%
  group_by(Subspecies) %>%
  summarise(TotalSites = sum(data_fold))

df$Proportional = df$data_fold/(TotalCountsFolded$TotalSites[match(df$Subspecies, TotalCountsFolded$Subspecies)]) #add column with proportion of sites

#Now plot
cbPalette = c("Generic" = "gray25", "Amur" = "#D55E00",  "Bengal" = "steelblue", "Malayan" = "#009E73", "Sumatran" = "gold3", "Indochinese" = "chocolate1", "Sanctuary" = "#CC79A7", "F1Wild" = "#867BCF", "Wild"="goldenrod2", "F2Wild"="coral2", "Captive"="cornflowerblue", "Generic"="darkmagenta", "Zoo"="darkseagreen3") #palette


propFoldedSFS_withGeneric = ggplot(df, aes(y=Proportional, x=bin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_fill_manual(values= cbPalette, name = "Sub-Species") + 
  scale_x_continuous(breaks = df$bin) + 
  labs(x= "SNP Frequency", y= "Proportion of SNPs") + 
  ggtitle("All Sub-Species (N=18)\nWhole Genome Folded SFS") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
        plot.title=element_text(size=24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18))

countsFoldedSFS_withGeneric = ggplot(df, aes(y=data_fold, x=bin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_fill_manual(values= cbPalette, name = "Sub-Species") + 
  scale_x_continuous(breaks = df$bin) + 
  labs(x= "SNP Frequency", y= "Count of SNPs") + 
  ggtitle("All Sub-Species (N=18)\nWhole Genome Folded SFS") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
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
        axis.text.y = element_text(size = 20),
        plot.title=element_text(size=24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18))



####unfolded sfs
setwd("~/TigerProject/sfs/Unfolded")
fnames = list.files(pattern = "\\N18.txt$")

unfoldedDF = rbindlist(sapply(fnames, read.delim, simplify = FALSE), use.names = TRUE, idcol = "Subspecies") %>%
  mutate(Subspecies = gsub("_.*$","",Subspecies))

#Plot Data
PlottingUnfolded = unfoldedDF %>%
  filter(FreqBin > 0 & FreqBin < 36) %>% #remove fixed stuff
  group_by(Subspecies, FreqBin) %>%
  summarise(WholeGenomeCounts = sum(sum)) %>%
  ungroup() 

TotalCounts = PlottingUnfolded %>%
  group_by(Subspecies) %>%
  summarise(TotalSites = sum(WholeGenomeCounts))

PlottingUnfolded$Proportional = PlottingUnfolded$WholeGenomeCounts/(TotalCounts$TotalSites[match(PlottingUnfolded$Subspecies, TotalCounts$Subspecies)])

#Plot only the bins starting with singletons
UnfoldedSFS = ggplot(PlottingUnfolded, aes(y=Proportional, x=FreqBin, fill=Subspecies)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))  +
  scale_x_continuous(breaks=1:35) + 
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
