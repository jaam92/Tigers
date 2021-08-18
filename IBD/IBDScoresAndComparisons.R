#Load libraries and data
library(tidyverse)
library(data.table)
library(ggpubr)

#Load ibdseq roh outputs and meta data
setwd("~/TigerProject/Rscripts/IBD")

popsDF = read_csv("~/TigerProject/individual_ids.csv")

#####IBD from truffle is reported in Mbp
df = read_delim("~/TigerProject/IBD/allChroms_truffle_allSubSpecies_calledPerSpecies.segments.coverage.bed", delim = "\t", col_names = c("chrom","start","end", "sample1","sample2", "typeIBD", "IBDLengthMb", "PropCovered"), col_types = "cnncccnn") %>%
  group_by(chrom,start,end,sample1,sample2,IBDLengthMb,typeIBD) %>%
  summarise_at(c("PropCovered"), sum, na.rm = TRUE) %>%
  ungroup()  

##Identify IBD greater than or equal to 2Mb and remove IBD where the avg. prop covered by SNPs is within 1 SD of the mean 
IBDgr2Mb = df[which(df$IBDLengthMb >= 2),][c("chrom","start","end", "sample1","sample2", "typeIBD", "IBDLengthMb", "PropCovered")]
z = data.table(IBDgr2Mb)
z[,ToKeep := abs(IBDgr2Mb$PropCovered - mean(IBDgr2Mb$PropCovered)) < sd(IBDgr2Mb$PropCovered)][ToKeep  == TRUE] #ID IBD to keep
dfIBD = subset(z, z$ToKeep == "TRUE") %>% #Subset out true IBD
  mutate(ToKeep = NULL,
         pop1 = popsDF$Subspecies[match(sample1, popsDF$Individual)],
         pop2 = popsDF$Subspecies[match(sample2, popsDF$Individual)]) 

groupScoreDF = dfIBD %>%
  filter(pop1 == pop2 & IBDLengthMb > 2 & IBDLengthMb <= 20) %>%
  group_by(pop1) %>%
  summarise(GroupScore = sum(as.numeric(IBDLengthMb))) %>% #sum up shared IBD for group
  ungroup() 

#count individuals per pop and get norm. constant
indivs = unique(c(unique(dfIBD$sample1), unique(dfIBD$sample2)))

IBDScoreDF = popsDF %>%
  filter(popsDF$Individual %in% indivs) %>%
  group_by(Subspecies) %>%
  count() %>%
  filter(n > 1) %>%
  mutate(normConstant = (choose(2*as.numeric(n), 2)) - as.numeric(n),
         GroupScore = groupScoreDF$GroupScore[match(Subspecies, groupScoreDF$pop1)],
         NormGroupScorePerMb = (GroupScore/normConstant)) %>%
  na.omit() #drop groups without ibd after filtering cut-offs

IBDScoreDF$RelToGeneric = IBDScoreDF$NormGroupScorePerMb/IBDScoreDF[grep("Unknown$", IBDScoreDF$Subspecies),]$NormGroupScorePerMb


write.table(IBDScoreDF, "~/TigerProject/IBD/IBDScores_includeRels.txt", quote = F, row.names = F, col.names = T, sep = "\t") #write out to file

##Plot data
cbPalette = c("Generic" = "gray25", "Generic* (Amur)" = "gray10" ,"Amur" = "#D55E00",  "Bengal" = "steelblue", "Malayan" = "#009E73", "Sumatran" = "gold3", "Indochinese" = "chocolate1", "Sanctuary" = "#CC79A7", "F1Wild" = "#867BCF", "Wild"="goldenrod2", "F2Wild"="coral2", "Captive"="cornflowerblue", "Unknown"="darkmagenta", "Zoo"="darkseagreen3", "South China" = "firebrick2") #palette

PlotIBDScores = ggplot(data = IBDScoreDF, aes(x=Subspecies, y=NormGroupScorePerMb)) +
  geom_point(aes(colour = Subspecies), size = 2) +
  labs(x="Sub-species", y="IBD Score (Mb)") +
  scale_colour_manual(name = "Sub-species", values = cbPalette) +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)) 

#plot relative to generic
mid = mean(IBDScoreDF$RelToGeneric)
PlotIBDRel2Generic = ggplot(IBDScoreDF, aes(y=RelToGeneric, x=Subspecies, fill=RelToGeneric)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  labs(y = "Within Population IBD Score (Mb) Relative to Generic\n Normalized by Sample Size", x="Sub-species") + 
  geom_hline(yintercept = mid, linetype="dashed", colour = "black") + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 10), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_gradient(name= "Fold-Enrichment", low="blue" ,high="red")

ggarrange(PlotIBDScores, PlotIBDRel2Generic, ncol = 2, align = 'h')
