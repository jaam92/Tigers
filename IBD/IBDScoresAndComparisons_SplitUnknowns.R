#Load libraries and data
library(tidyverse)
library(data.table)
library(ggpubr)

#Load ibdseq roh outputs and meta data
setwd("~/TigerProject/IBD")

popsDF = read_csv("~/TigerProject/IndivFiles/individual_ids.csv") %>%
  mutate(combo = ifelse(Subspecies == "Generic", paste(Subspecies, "-", Phenotype, sep = ""), Subspecies))


uSub2 = read.csv("~/TigerProject/IndivFiles/Unrelateds_basedOnSubspecies2_highCov.csv") #unrelated based on subspecies 2, split Generics

#####IBD from truffle is reported in Mbp
df = read_delim("~/TigerProject/IBD/allChroms_truffle_allSubSpecies_calledPerSpecies.segments.coverage.bed", delim = "\t", col_names = c("chrom","start","end", "sample1","sample2", "typeIBD", "IBDLengthMb", "PropCovered"), col_types = "cnncccnn") %>%
  group_by(chrom,start,end,sample1,sample2,IBDLengthMb,typeIBD) %>%
  summarise_at(c("PropCovered"), sum, na.rm = TRUE) %>%
  ungroup()  

##Identify IBD greater than or equal to 2Mb and remove IBD where the avg. prop covered by SNPs is within 1 SD of the mean 
IBDgr2Mb = df[which(df$IBDLengthMb >= 2),][c("chrom","start","end", "sample1","sample2", "typeIBD", "IBDLengthMb", "PropCovered")]
z = data.table(IBDgr2Mb)
z[,ToKeep := abs(IBDgr2Mb$PropCovered - mean(IBDgr2Mb$PropCovered)) < sd(IBDgr2Mb$PropCovered)][ToKeep  == TRUE] #ID IBD to keep


#split on Subspecies 2
dfIBD = z %>%
  filter(ToKeep == "TRUE" & sample1%in%uSub2$Individual & sample2%in%uSub2$Individual) %>%
  mutate(ToKeep = NULL,
         pop1 = popsDF$combo[match(sample1, popsDF$Individual)],
         pop2 = popsDF$combo[match(sample2, popsDF$Individual)])

groupScoreDF= dfIBD %>%
  filter(pop1 == pop2 & IBDLengthMb > 2 & IBDLengthMb <= 20) %>%
  group_by(pop1) %>%
  summarise(GroupScore = sum(as.numeric(IBDLengthMb))) %>% #sum up shared IBD for group
  ungroup() 

#count individuals per pop and get norm. constant
indivs = unique(c(unique(dfIBD$sample1), unique(dfIBD$sample2)))

IBDScoreDF = popsDF %>%
  filter(Individual %in% indivs) %>%
  group_by(combo) %>%
  count() %>%
  filter(n > 1) %>%
  mutate(normConstant = (choose(2*as.numeric(n), 2)) - as.numeric(n),
         GroupScore = groupScoreDF$GroupScore[match(combo, groupScoreDF$pop1)],
         NormGroupScorePerMb = (GroupScore/normConstant)) %>%
  na.omit() #drop groups without ibd after filtering cut-offs

IBDScoreDF$RelToGeneric = IBDScoreDF$NormGroupScorePerMb/IBDScoreDF[grep("Generic-Orange$", IBDScoreDF$combo),]$NormGroupScorePerMb

write.table(IBDScoreDF, "~/TigerProject/IBD/IBDScores_subspecies2Unrelateds.txt", quote = F, row.names = F, col.names = T, sep = "\t") #write out to file


##Plot data
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Sumatran" = "cornflowerblue", "Indochinese" = "gold4", "South China" = "plum", "Generic"="gray25", "Generic-Orange" = "#CC79A7", "Generic-SnowWhite" = "#867BCF", "Generic-Golden"="darkseagreen3", "Generic-White"="cornflowerblue")#palette

PlotIBDScores = ggplot(data = IBDScoreDF, aes(x=combo, y=NormGroupScorePerMb)) +
  geom_point(aes(colour = combo), size = 2) +
  labs(x="Sub-species", y="IBD Score (Mb)") +
  scale_colour_manual(name = "Sub-species", values = cbPalette) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 18), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18)) 

#plot relative to generic
mid = mean(IBDScoreDF$RelToGeneric)
PlotIBDRel2Generic = ggplot(IBDScoreDF, aes(y=RelToGeneric, x=combo, fill=RelToGeneric)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  labs(y = "Within Population IBD Score (Mb) Relative to Generic-Orange\n Normalized by Sample Size", x="Sub-species") + 
  geom_hline(yintercept = mid, linetype="dashed", colour = "black") + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 18), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_gradient(name= "Fold-Enrichment", low="blue" ,high="red")

ggarrange(PlotIBDScores, PlotIBDRel2Generic, ncol = 2, align = 'h')
