#Load libraries and data
library(tidyverse)
library(data.table)
library(ggpubr)

#Load ibdseq roh outputs and meta data
setwd("~/Documents/Tigers/IBD")

popsDF = read_delim("~/Documents/Tigers/IndivFiles/TableX-SampleDetails.txt") 
unrelateds = read_delim("~/Documents/Tigers/IndivFiles/N10-N6_unrelateds.txt", delim = "\t")

# ###Only needed to run this to generate the files initially
# df = read_delim("~/Documents/Tigers/IBD/allChroms_truffle_allSubSpecies_calledPerSpecies.segments.coverage.bed", delim = "\t", col_names = c("chrom","start","end", "sample1","sample2", "typeIBD", "IBDLengthMb", "PropCovered"), col_types = "cnncccnn") %>%
#   group_by(chrom,start,end,sample1,sample2,IBDLengthMb,typeIBD) %>%
#   summarise_at(c("PropCovered"), sum, na.rm = TRUE) %>%
#   ungroup() 
# 
# ##Identify IBD greater than or equal to 2Mb and remove IBD where the avg. prop covered by SNPs is within 1 SD of the mean
# IBDgr2Mb = df[which(df$IBDLengthMb >= 2),][c("chrom","start","end", "sample1","sample2", "typeIBD", "IBDLengthMb", "PropCovered")]
# z = data.table(IBDgr2Mb)
# z[,ToKeep := abs(IBDgr2Mb$PropCovered - mean(IBDgr2Mb$PropCovered)) < sd(IBDgr2Mb$PropCovered)][ToKeep  == TRUE] #ID IBD to keep

#write.table(z, "~/Documents/Tigers/IBD/allIBD_propCoveredwithin1SDMean_gr2Mb_allChroms_highCov_runSpeciesSep_truffle.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#write.table(z %>% filter(ToKeep == "TRUE"), "~/Documents/Tigers/IBD/TrueIBD_propCoveredwithin1SDMean_gr2Mb_allChroms_highCov_runSpeciesSep_truffle.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#IBD segments to keep (only include unrelateds)
z = read_delim("~/Documents/Tigers/IBD/allIBD_propCoveredwithin1SDMean_gr2Mb_allChroms_highCov_runSpeciesSep_truffle.txt", delim = "\t") 

dfIBD = z %>%
  filter(ToKeep == "TRUE" & sample1 %in% unrelateds$Sample & sample2 %in% unrelateds$Sample) %>%
  mutate(ToKeep = NULL,
         pop1 = popsDF$Subspecies_GroupID_Corrected[match(sample1, popsDF$Sample)],
         pop2 = popsDF$Subspecies_GroupID_Corrected[match(sample2, popsDF$Sample)]) 

groupScoreDF = dfIBD %>%
  filter(pop1 == pop2 & IBDLengthMb > 2 & IBDLengthMb <= 20) %>%
  group_by(pop1) %>%
  summarise(GroupScore = sum(as.numeric(IBDLengthMb))) %>% #sum up shared IBD for group
  ungroup() 


#count individuals per pop and get norm. constant
indivs = unique(c(unique(dfIBD$sample1), unique(dfIBD$sample2)))

IBDScoreDF = popsDF %>%
  filter(Sample %in% indivs) %>%
  group_by(Subspecies_GroupID_Corrected) %>%
  count() %>%
  filter(n > 1) %>%
  mutate(normConstant = (choose(2*as.numeric(n), 2)) - as.numeric(n),
         GroupScore = groupScoreDF$GroupScore[match(Subspecies_GroupID_Corrected, groupScoreDF$pop1)],
         NormGroupScorePerMb = (GroupScore/normConstant),
         Subspecies2 = factor(Subspecies_GroupID_Corrected, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran'))) %>%
  na.omit() #drop groups without ibd after filtering cut-offs

IBDScoreDF$RelToGeneric = IBDScoreDF$NormGroupScorePerMb/IBDScoreDF[grep("Generic$", IBDScoreDF$Subspecies2),]$NormGroupScorePerMb
#write.table(IBDScoreDF, "~/Documents/Tigers/IBD/IBDScores_subspecies2_unrelatedSamps.txt", quote = F, row.names = F, col.names = T, sep = "\t") #write out to file


sharing = dfIBD %>%
  mutate(key = paste(pmin(sample1, sample2), pmax(sample1, sample2), sep = "-")) %>%
  group_by(pop1, key) %>%
  summarise(GroupScore = sum(as.numeric(IBDLengthMb))) %>% #sum up shared ROH within IBD for group
  ungroup() %>%
  mutate(pop1 = factor(pop1, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran')))

##Plot data
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Indochinese" = "gold4", "South China" = "plum", "Sumatran" = "cornflowerblue", "Generic"="gray25")#palette
cbPalette_expanded = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Indochinese" = "gold4", "South China" = "plum", "Sumatran" = "cornflowerblue", "Generic"="gray25", "Generic-Orange" = "#CC79A7", "Generic-SnowWhite" = "#867BCF", "Generic-Golden"="darkseagreen3", "Generic-White"="cornflowerblue")#palette

PlotIBDScores = ggplot(data = IBDScoreDF, aes(x=Subspecies2, y=NormGroupScorePerMb)) +
  geom_point(aes(colour = Subspecies2), size = 2) +
  labs(x="Subspecies", y="IBD Score (Mb)") +
  scale_colour_manual(name = "Subspecies", values = cbPalette) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 18), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18)) 

#plot relative to generic
mid = mean(IBDScoreDF$RelToGeneric)
PlotIBDRel2Generic = ggplot(IBDScoreDF, aes(y=RelToGeneric, x=Subspecies2, fill=RelToGeneric)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  labs(y = "Within Population IBD Score (Mb) Relative to Generic\n Normalized by Sample Size", x="Subspecies") + 
  geom_hline(yintercept = mid, linetype="dashed", colour = "black") + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 18), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_gradient(name= "Fold-Enrichment", low="blue" ,high="red")

print(PlotIBDScores)
ggarrange(PlotIBDScores + theme(legend.position = "none"), PlotIBDRel2Generic, ncol = 2, align = 'h')

pairwiseSharing = ggplot(data = sharing, aes(x=pop1, y=GroupScore, fill=pop1)) +
  geom_jitter(height = 0, width = 0.1, shape=23, size = 3) +
  labs(x="Subspecies", y="Pairwise sharing\n between individuals (Mb)") +
  scale_fill_manual(name = "Subspecies", values = cbPalette) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 18), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18)) 

print(pairwiseSharing)
