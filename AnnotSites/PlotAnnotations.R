#Load Libraries
library(tidyverse)
library(ggpubr)

#Plotting fxn and color palette
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Sumatran" = "cornflowerblue", "Indochinese" = "gold4", "South China" = "plum", "Generic"="gray25")#palette

plotFxn = function(dataFrame, x_axisCol, y_axisCol, y_axisTitle) {
  RaincloudWithBoxPlot = ggplot(dataFrame, aes(x_axisCol, y_axisCol, fill=Subspecies)) + 
    ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) + 
    geom_boxplot(width = .2, outlier.shape = NA) + 
    geom_jitter(width = .05, alpha = .5) +
    coord_flip() +
    guides(fill = "none", colour = "none") +
    scale_fill_manual(values = cbPalette) + 
    labs(x="Subspecies",y=paste0(y_axisTitle)) + 
    theme_bw() + 
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16), 
          axis.text.y = element_text(size = 16), 
          plot.title = element_text(size = 16, hjust = 0.5), 
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 14))
  return(RaincloudWithBoxPlot)
}

#Read file in 
removed = c('T18', 'T5', 'T10', 'SRR5591010', 'SRR5612311', 'SRR5612312')
removed_plusoutlier = c('T18', 'T5', 'T10', 'SRR5591010', 'SRR5612311', 'SRR5612312', 'BEN_NE2', 'GEN1') #last two individuals are outliers on the plot
PlotDF = read_delim("~/Documents/Tigers/AnnotSites/GTAnnotationCountResults_Nov2022_Tigers_addSIFT.txt", delim = "\t") %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  mutate(CallableSites = LineCount - Missing,
         Subspecies = factor(Subspecies, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran'))) %>%
  filter(!ID %in% removed_plusoutlier & Missing < 2500)
  #filter(!ID %in% removed & Missing < 2500) 
  


####Make Everything proportional 
PropPlotDF = PlotDF[,c(1:11, 20:64, 71:106)] %>%
  mutate_at(vars(LOF_CountAlleles:PutNeu_CountVariants_ROH_TypeC), list(~./PlotDF$CallableSites)) #Just run replace on plotting fxns with regular df to make everything proportional

####Make Scaled Count Alleles
scaleCalls = mean(PlotDF$CallableSites)
ScaledPlotDF = PropPlotDF %>%
  mutate_at(vars(LOF_CountAlleles:PutNeu_CountVariants_ROH_TypeC), list(~.*scaleCalls)) %>%
  mutate_if(is.numeric, ceiling) 

#write.table(ScaledPlotDF, file = "~/Documents/Tigers/AnnotSites/scaledGTAnnotationCountResults_Nov2022_Tigers.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

###Plot Supplementary Figures putative neutral and deleterious
CountDerHom_SY = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies ,y_axisCol = ScaledPlotDF$SY_CountDerHom, y_axisTitle = "Count derived neutral homozygotes")
CountDerHom_NS = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies ,y_axisCol = ScaledPlotDF$NS_CountDerHom, y_axisTitle = "Count derived deleterious homozygotes")

CountVar_SY = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies ,y_axisCol = ScaledPlotDF$SY_CountVariants, y_axisTitle = "Count derived neutral variants")
CountVar_NS = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies ,y_axisCol = ScaledPlotDF$NS_CountVariants, y_axisTitle = "Count derived deleterious variants")

CountAllele_SY = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies ,y_axisCol = ScaledPlotDF$SY_CountAlleles, y_axisTitle = "Count derived neutral alleles")
CountAllele_NS = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies ,y_axisCol = ScaledPlotDF$NS_CountAlleles, y_axisTitle = "Count derived deleterious alleles")


#Neutral
SY = ggarrange( CountDerHom_SY + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count Homozygotes"),
                    CountVar_SY + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count Variants"),
                    CountAllele_SY + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count Alleles"),
                    align = 'hv',
                    #labels = c("A", "B", "C"),
                    hjust = -1,
                    nrow = 1)
SYAnnot = annotate_figure(SY, left = text_grob("Synonymous", size=16, rot = 90, hjust = 0.5))

#Deleterious
NS = ggarrange( CountDerHom_NS + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                    CountVar_NS + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                    CountAllele_NS + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                    align = 'hv',
                    #labels = c("D", "E", "F"),
                    hjust = -1,
                    nrow = 1)
NSAnnot = annotate_figure(NS, left = text_grob("Nonsynonymous", size=16, rot = 90, hjust = 0.5))

#plot deleterious and neutral
ggarrange(SYAnnot, NSAnnot, nrow = 2)

#arrange the three scaled derived allele count plots in a single row
scaledCountDerHom_PutNeu = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies ,y_axisCol = ScaledPlotDF$PutNeuSIFT_CountDerHom, y_axisTitle = "Count derived neutral homozygotes")
scaledCountDerHom_PutDel = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies ,y_axisCol = ScaledPlotDF$PutDelSIFT_CountDerHom, y_axisTitle = "Count derived deleterious homozygotes")

scaledCountVar_PutNeu = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies ,y_axisCol = ScaledPlotDF$PutNeuSIFT_CountVariants, y_axisTitle = "Count derived neutral variants")
scaledCountVar_PutDel = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies ,y_axisCol = ScaledPlotDF$PutDelSIFT_CountVariants, y_axisTitle = "Count derived deleterious variants")

scaledCountAllele_PutNeu = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies ,y_axisCol = ScaledPlotDF$PutNeuSIFT_CountAlleles, y_axisTitle = "Count derived neutral alleles")
scaledCountAllele_PutDel = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies ,y_axisCol = ScaledPlotDF$PutDelSIFT_CountAlleles, y_axisTitle = "Count derived deleterious alleles")

#Neutral
scaledPutNeu = ggarrange(scaledCountDerHom_PutNeu + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count Homozygotes"),
                          scaledCountVar_PutNeu + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count of Variants"),
                          scaledCountAllele_PutNeu + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count of Alleles"),
                    align = 'hv',
                    #labels = c("A", "B", "C"),
                    hjust = -1,
                    nrow = 1)
scaledPutNeuAnnot = annotate_figure(scaledPutNeu, left = text_grob("Neutral", size=16,rot = 90, hjust = 0.5))

#Deleterious
scaledPutDel = ggarrange( scaledCountDerHom_PutDel + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                          scaledCountVar_PutDel + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                          scaledCountAllele_PutDel + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                    align = 'hv',
                    #labels = c("D", "E", "F"),
                    hjust = -1,
                    nrow = 1)
scaledPutDelAnnot = annotate_figure(scaledPutDel, left = text_grob("Deleterious", size=16,rot = 90, hjust = 0.5))

#plot deleterious and neutral
ggarrange(scaledPutNeuAnnot, scaledPutDelAnnot, nrow = 2)

#Compute pvalues
pairwise.wilcox.test(ScaledPlotDF$PutDelSIFT_CountDerHom, ScaledPlotDF$Subspecies, p.adj = "bonf")$p.value
pairwise.wilcox.test(ScaledPlotDF$PutDelSIFT_CountVariants, ScaledPlotDF$Subspecies, p.adj = "bonf")$p.value
pairwise.wilcox.test(ScaledPlotDF$PutDelSIFT_CountAlleles, ScaledPlotDF$Subspecies, p.adj = "bonf")$p.value




######making Odds-ratios
test = ScaledPlotDF %>%
  group_by(Subspecies) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE) %>%
  mutate(PutDel_ROH = PutDel_CountDerHom_ROH_TypeA + PutDel_CountDerHom_ROH_TypeB + PutDel_CountDerHom_ROH_TypeC,
         PutNeu_ROH = PutNeu_CountDerHom_ROH_TypeA + PutNeu_CountDerHom_ROH_TypeB + PutNeu_CountDerHom_ROH_TypeC,
         PutDel_nonROH = PutDel_CountDerHom_nonROH_TypeA + PutDel_CountDerHom_nonROH_TypeB + PutDel_CountDerHom_nonROH_TypeC,
         PutNeu_nonROH = PutNeu_CountDerHom_nonROH_TypeA + PutNeu_CountDerHom_nonROH_TypeB + PutNeu_CountDerHom_nonROH_TypeC,
         NS_ROH = NS_CountDerHom_ROH_TypeA + NS_CountDerHom_ROH_TypeB + NS_CountDerHom_ROH_TypeC,
         SY_ROH = SY_CountDerHom_ROH_TypeA + SY_CountDerHom_ROH_TypeB,
         NS_nonROH = NS_CountDerHom_nonROH_TypeA + NS_CountDerHom_nonROH_TypeB + NS_CountDerHom_nonROH_TypeC,
         SY_nonROH = SY_CountDerHom_nonROH_TypeA + SY_CountDerHom_nonROH_TypeB ,
         OR_DerHom_ROHvnonROH_PutDelPutNeu = (PutNeu_nonROH*PutDel_ROH)/(PutNeu_ROH*PutDel_nonROH),
         OR_DerHom_ROHvnonROH_SYNS = (SY_nonROH*NS_ROH)/(SY_ROH*NS_nonROH),
         PutDel_rel_PutNeu_DerHom = PutDelSIFT_CountDerHom/PutNeuSIFT_CountDerHom,
         NS_rel_SY_DerHom = NS_CountDerHom/SY_CountDerHom)


#subset data
x = test %>% 
  select(Subspecies, OR_DerHom_ROHvnonROH_PutDelPutNeu) %>% 
  mutate(Type = "Putatively Deleterious Relative to Putatively Neutral") %>%
  rename("Analysis" = "OR_DerHom_ROHvnonROH_PutDelPutNeu")

y = test %>% 
  select(Subspecies, OR_DerHom_ROHvnonROH_SYNS) %>% 
  mutate(Type = "Nonsynonymous Relative to Synonymous") %>% 
  rename("Analysis" = "OR_DerHom_ROHvnonROH_SYNS")

z = rbind.data.frame(x,y)


j = data.frame() 
for (i in unique(test$Subspecies)) {
  
  f = test %>% 
    filter(Subspecies == i)
  
  l = matrix(c(as.numeric(f["SY_nonROH"]), as.numeric(f["SY_ROH"]), as.numeric(f["NS_nonROH"]), as.numeric(f["NS_ROH"])),
             nrow=2,
             byrow=TRUE)
  
  m = cbind.data.frame(f["OR_DerHom_ROHvnonROH_SYNS"], 
                       fisher.test(l)$conf.int[1],
                       fisher.test(l)$conf.int[2],
                       fisher.test(l)$p.value,
                       f["Subspecies"])
  
  colnames(m) = c("OR", "left", "right","pvalue", "Subspecies")
  
  m$Analysis = "Nonsynonymous Relative to Synonymous"
  m$significant = ifelse(m$pvalue <= 0.05, "yes","no") 
  
  j = rbind.data.frame(m, j)
  
}

#Putatively del and neu
for (i in unique(test$Subspecies)) {
  
  f = test %>% 
    filter(Subspecies == i)
  
  l = matrix(c(as.numeric(f["PutNeu_nonROH"]), as.numeric(f["PutNeu_ROH"]), as.numeric(f["PutDel_nonROH"]), as.numeric(f["PutDel_ROH"])),
             nrow=2,
             byrow=TRUE)
  
  m = cbind.data.frame(f["OR_DerHom_ROHvnonROH_PutDelPutNeu"], 
                       fisher.test(l)$conf.int[1],
                       fisher.test(l)$conf.int[2],
                       fisher.test(l)$p.value,
                       f["Subspecies"])
  
  colnames(m) = c("OR", "left", "right","pvalue", "Subspecies")
  
  m$Analysis = "Putatively Deleterious Relative to Putatively Neutral"
  m$significant = ifelse(m$pvalue <= 0.05, "yes","no") 
  
  j = rbind.data.frame(m, j)
  
}


#plot
ggplot(j, aes(x=Subspecies, y=OR, color=significant)) +
  facet_grid(~Analysis) +
  geom_errorbar(aes(ymin=left, ymax=right), colour="gray40", width=.2) + 
  geom_point() + 
  coord_flip() +  
  scale_colour_manual(values = c("yes"= "red", "no"="black"), guide="none") +
  labs(x="Subspecies", y="Odds-ratio") + 
  ylim(0,3.5) +
  geom_hline(yintercept = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 24), 
        axis.text.y = element_text(size = 24), 
        plot.title = element_text(size = 24, hjust = 0.5), 
        axis.title = element_text(size = 24),
        strip.text = element_text(size = 24),
        legend.position = "none")
