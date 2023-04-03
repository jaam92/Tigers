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
PlotDF = read_delim("~/Documents/Tigers/AnnotSites/GTAnnotationCountResults_Nov2022_Tigers.txt", delim = "\t") %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  mutate(CallableSites = LineCount - Missing,
         Subspecies = factor(Subspecies, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran'))) %>%
  filter(!ID %in% removed_plusoutlier & Missing < 2500)
  #filter(!ID %in% removed & Missing < 2500) 
  


####Make Everything proportional 
PropPlotDF = PlotDF[,c(1:11, 20:25)] %>%
  mutate_at(vars(LOF_CountAlleles:SY_CountVariants), list(~./PlotDF$CallableSites)) #Just run replace on plotting fxns with regular df to make everything proportional

####Make Scaled Count Alleles
scaleCalls = mean(PlotDF$CallableSites)
ScaledPlotDF = PropPlotDF %>%
  mutate_at(vars(LOF_CountAlleles:SY_CountVariants), list(~.*scaleCalls)) %>%
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




####combine with FSNP and FROH
source(file = "~/Documents/Tigers/ROH/plotROH_garlic.R")

mergedPlotDF = ScaledPlotDF %>%
  mutate(FROH = x$FROH[match(ID, x$INDV)],
         FSNP = all_nodups_full_highCov$FSNP[match(ID, all_nodups_full_highCov$Sample)],
         Heterozygosity = all_nodups_full_highCov$Heterozygosity[match(ID, all_nodups_full_highCov$Sample)]) %>%
  na.omit(FSNP) %>%
  na.omit(FROH) %>%
  na.omit(Heterozygosity)


putDelAllele_FROH = ggplot(mergedPlotDF, aes(x=FROH, y=PutDelSIFT_CountAlleles)) +
  geom_point(aes(colour = Subspecies), size = 2) +
  geom_smooth(method = "lm") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), color = "red", geom = "label") + #this is not the adjusted r2
  facet_wrap(.~Subspecies) +
  scale_colour_manual(name = "Subspecies", values = cbPalette) +
  labs(x=expression(F[ROH]), y="Count derived deleterious alleles") +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14))

putDelDerHom_FROH = ggplot(mergedPlotDF, aes(x=FROH, y=PutDelSIFT_CountAlleles)) +
  geom_point(aes(colour = Subspecies), size = 2) +
  geom_smooth(method = "lm") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), color = "red", geom = "label") + #this is not the adjusted r2
  facet_wrap(.~Subspecies) +
  scale_colour_manual(name = "Subspecies", values = cbPalette) +
  labs(x=expression(F[ROH]), y="Count derived deleterious homozygotes") +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)) 

putDelDerHom_FSNP = ggplot(mergedPlotDF, aes(x=FSNP, y=PutDelSIFT_CountAlleles)) +
  geom_point(aes(colour = Subspecies), size = 2) +
  geom_smooth(method = "lm") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), color = "red", geom = "label") + #this is not the adjusted r2
  facet_wrap(.~Subspecies) +
  scale_colour_manual(name = "Subspecies", values = cbPalette) +
  labs(x=expression(F[SNP]), y="Count derived deleterious homozygotes") +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)) 

putDelAllele_FSNP = ggplot(mergedPlotDF, aes(x=FSNP, y=PutDelSIFT_CountAlleles)) +
  geom_point(aes(colour = Subspecies), size = 2) +
  geom_smooth(method = "lm") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), color = "red", geom = "label") + #this is not the adjusted r2
  facet_wrap(.~Subspecies) +
  scale_colour_manual(name = "Subspecies", values = cbPalette) +
  labs(x=expression(F[SNP]), y="Count derived deleterious alleles") +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14))

ggarrange(putDelAllele_FROH, putDelDerHom_FROH, common.legend = TRUE)
ggarrange(putDelAllele_FSNP, putDelDerHom_FSNP, common.legend = TRUE)


