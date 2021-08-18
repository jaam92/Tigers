#Load libraries and data
library(tidyverse)
library(ggpubr)

#Load bcftools roh outputs and meta data
setwd("~/Documents/TigerProject/Rscripts/")
highCoverage = read_delim("~/Documents/TigerProject/IndivFiles/highcov-info.txt", delim = "\t") %>%
  select(-c(Other_information)) %>%
  arrange(Species)

df = read.delim("TrueROH_propCoveredwithin1SDMean_gr100kb_gr50SNPs_allChroms_highCov_runSpeciesSep.txt")  %>%
  rename("Sample"="INDV") %>%
  left_join(highCoverage, by = "Sample") %>%
  mutate(Origin = ifelse(Origin == "Sanctuary" | Origin == "Zoo", "Captive", Origin),
         Origin = ifelse(Origin == "F1Wild" | Origin == "F2Wild", "Wild", Origin))

#Plottting DF
plotGeneric = df %>%
  group_by(Sample, Species, Origin) %>%
  #filter(AUTO_LEN > 10^6) %>%
  summarise_at(c("AUTO_LEN"), sum, na.rm=TRUE) %>%
  ungroup() 

cbPalette = c("Generic" = "gray25", "Amur" = "#D55E00",  "Bengal" = "steelblue", "Malayan" = "#009E73", "Sumatran" = "gold3", "Indochinese" = "chocolate1", "Sanctuary" = "#CC79A7", "F1Wild" = "#867BCF", "Wild"="goldenrod2", "F2Wild"="coral2", "Captive"="cornflowerblue", "Unknown"="darkmagenta", "Zoo"="darkseagreen3") #palette

ROH_Species = ggplot(plotGeneric , aes(y=Species, x=AUTO_LEN/10^6)) +
  geom_boxplot(outlier.shape = NA) + #remove outlier points and only use jitter
  geom_jitter(aes(colour = Species, shape = Origin), size = 2) +
  labs(y="Species", x="Length of Genome in an ROH (Mb) per Individual") +
  scale_colour_manual(name = "Species", values = cbPalette) + 
  guides(colour = guide_legend(title = "Species", order = 1), 
         shape = guide_legend(title = "Origin", order = 2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        plot.title=element_text(size=24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))


###Ranges
#Generate ROH Ranges
#Split into length intervals, find length for each indiv, add column in Gb, add range column
aggIntROH = df %>% 
  select(AUTO_LEN, Sample) %>% 
  filter(AUTO_LEN >= 100000 & AUTO_LEN < 1000000) %>% 
  group_by(Sample) %>% 
  summarise(totalLen = sum(AUTO_LEN)) %>% 
  mutate(totalLenGb = totalLen/10^9, 
         Range = "[0.1-1)Mb") %>% 
  as.data.frame()

aggLongROH = df %>% 
  select(AUTO_LEN, Sample) %>% 
  filter(AUTO_LEN >= 1000000 & AUTO_LEN < 10000000) %>% 
  group_by(Sample) %>% 
  summarise(totalLen = sum(AUTO_LEN)) %>% 
  mutate(totalLenGb = totalLen/10^9, 
         Range = "[1-10)Mb")  %>% 
  as.data.frame()

aggVLongROH = df %>% 
  select(AUTO_LEN, Sample) %>% 
  filter(AUTO_LEN >= 10000000) %>% 
  group_by(Sample) %>% 
  summarise(totalLen = sum(AUTO_LEN)) %>% 
  mutate(totalLenGb = totalLen/10^9, 
         Range = "[10-63)Mb") %>% 
  as.data.frame()

#Merge interval data frames
FinalDF = rbind(aggIntROH, aggLongROH, aggVLongROH)
FinalDF$Range = factor(FinalDF$Range, levels=c("[10-63)Mb","[1-10)Mb","[0.1-1)Mb"))
FinalDF$Species = highCoverage$Species[match(FinalDF$Sample, highCoverage$Sample)]
FinalDF$Sample = factor(FinalDF$Sample, levels = highCoverage$Sample)

#Function to plot ROH ranges with all ROH greater than 100Kb
plotROHRanges = function(dataFrame){
  ggplot(dataFrame, aes(x=Sample, y=totalLen/10^9, fill=Range)) + 
    geom_bar(stat="identity") + 
    theme_bw() + 
    coord_flip() + 
    scale_fill_manual(values = c("[0.1-1)Mb"= "bisque3", 
                                 "[1-10)Mb" = "darkgoldenrod",
                                 "[10-63)Mb" = "indianred4"), 
                      breaks = c("[0.1-1)Mb",
                                 "[1-10)Mb", 
                                 "[10-63)Mb"), 
                      name = "Range") + 
    labs(x ="Individual", y = "ROH length per bin (Gb)") + 
    theme(axis.text.x = element_text(size = 20), 
          axis.text.y = element_text(size = 10),
          plot.title=element_text(size=24, face = "bold", hjust=0.5), 
          axis.title=element_text(size=20, face = "bold"),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=18))
}


#ROH Ranges plots
ROHRangePlot = plotROHRanges(FinalDF)
