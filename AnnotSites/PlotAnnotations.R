#Load Libraries
source("~/DogProject_Jaz/LocalRscripts/OMIA/R_rainclouds.R")
library(tidyverse)
library(ggpubr)

#Plotting fxn and color palette
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Sumatran" = "cornflowerblue", "Indochinese" = "gold4", "South China" = "plum", "Unknown"="gray25", "Unknown-Orange" = "#CC79A7", "Unknown-SnowWhite" = "#867BCF", "Unknown-Golden"="darkseagreen3", "Unknown-White"="cornflowerblue")#palette #palette

plotFxn = function(dataFrame, x_axisCol, y_axisCol, y_axisTitle) {
  RaincloudWithBoxPlot = ggplot(dataFrame, aes(x=x_axisCol, y=y_axisCol, colour=x_axisCol)) +
    geom_flat_violin(size=1, position = position_nudge(x = .25, y = 0),adjust =2, trim = FALSE) +
    geom_point(aes(x = x_axisCol, y = y_axisCol, 
                   colour = x_axisCol), position = position_jitter(width = .05),
               size = 1, shape = 20) +
    geom_boxplot(aes(x = x_axisCol, y = y_axisCol, fill = x_axisCol),
                 outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
    coord_flip() +
    guides(fill = FALSE, colour = FALSE) +
    scale_colour_manual(values = cbPalette) +
    scale_fill_manual(values = cbPalette) + 
    labs(x="Subspecies",y=paste0(y_axisTitle)) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size=20, angle = 45, vjust=0.7), 
          axis.text.y = element_text(size=20), 
          plot.title=element_text(size=24, hjust = 0.5, face = "bold"), 
          axis.title=element_text(size=24),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=20))
  return(RaincloudWithBoxPlot)
}

#Read file in 
PlotDF = read_delim("~/TigerProject/AnnotSites/GTAnnotationCountResults_Sept2021_Tigers.txt", delim = "\t") %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  mutate(CallableSites = LineCount - Missing)


####Make Everything proportional 
PropPlotDF = PlotDF[,c(1:12, 21:26)] %>%
  mutate_at(vars(LOF_CountAlleles: PutNeuSIFT_CountVariants), list(~./PlotDF$CallableSites)) #Just run replace on plotting fxns with regular df to make everything proportional

####Make Scaled Count Alleles
scaleCalls = mean(PlotDF$CallableSites)
ScaledPlotDF = PropPlotDF %>%
  mutate_at(vars(LOF_CountAlleles: PutNeuSIFT_CountVariants), list(~.*scaleCalls))
#write.table(ScaledPlotDF, file = "~/TigerProject/AnnotSites/scaledGTAnnotationCountResults_Sept2021_Tigers.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

###Plot Supplementary Figures putativeley neutral and deleterious
CountDerHom_PutNeu = plotFxn(dataFrame = PlotDF, x_axisCol = PlotDF$Subspecies ,y_axisCol = PlotDF$PutNeuSIFT_CountDerHom, y_axisTitle = "Count derived neutral homozygotes")
CountDerHom_PutDel = plotFxn(dataFrame = PlotDF, x_axisCol = PlotDF$Subspecies ,y_axisCol = PlotDF$PutDelSIFT_CountDerHom, y_axisTitle = "Count derived deleterious homozygotes")

CountVar_PutNeu = plotFxn(dataFrame = PlotDF, x_axisCol = PlotDF$Subspecies ,y_axisCol = PlotDF$PutNeuSIFT_CountVariants, y_axisTitle = "Count derived neutral variants")
CountVar_PutDel = plotFxn(dataFrame = PlotDF, x_axisCol = PlotDF$Subspecies ,y_axisCol = PlotDF$PutDelSIFT_CountVariants, y_axisTitle = "Count derived deleterious variants")

CountAllele_PutNeu = plotFxn(dataFrame = PlotDF, x_axisCol = PlotDF$Subspecies ,y_axisCol = PlotDF$PutNeuSIFT_CountAlleles, y_axisTitle = "Count derived neutral alleles")
CountAllele_PutDel = plotFxn(dataFrame = PlotDF, x_axisCol = PlotDF$Subspecies ,y_axisCol = PlotDF$PutDelSIFT_CountAlleles, y_axisTitle = "Count derived deleterious alleles")


#Neutral
putNeu = ggarrange( CountDerHom_PutNeu + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count Homozygotes"),
                    CountVar_PutNeu + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count Variants"),
                    CountAllele_PutNeu + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count Alleles"),
                    align = 'hv',
                    labels = c("A", "B", "C"),
                    hjust = -1,
                    nrow = 1)
putNeuAnnot = annotate_figure(putNeu, left = text_grob("Neutral", size=24, face="bold",rot = 90, hjust = 0.5))

#Deleterious
putDel = ggarrange( CountDerHom_PutDel + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                    CountVar_PutDel + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                    CountAllele_PutDel + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                    align = 'hv',
                    labels = c("D", "E", "F"),
                    hjust = -1,
                    nrow = 1)
putDelAnnot = annotate_figure(putDel, left = text_grob("Deleterious", size=24, face="bold",rot = 90, hjust = 0.5))

#plot deleterious and neutral
ggarrange(putNeuAnnot, putDelAnnot, nrow = 2)

#arrange the three scaled derived allele count plots in a single row
scaledCountDerHom_PutNeu = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies2 ,y_axisCol = ScaledPlotDF$PutNeuSIFT_CountDerHom, y_axisTitle = "Count derived neutral homozygotes")
scaledCountDerHom_PutDel = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies2 ,y_axisCol = ScaledPlotDF$PutDelSIFT_CountDerHom, y_axisTitle = "Count derived deleterious homozygotes")

scaledCountVar_PutNeu = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies2 ,y_axisCol = ScaledPlotDF$PutNeuSIFT_CountVariants, y_axisTitle = "Count derived neutral variants")
scaledCountVar_PutDel = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies2 ,y_axisCol = ScaledPlotDF$PutDelSIFT_CountVariants, y_axisTitle = "Count derived deleterious variants")

scaledCountAllele_PutNeu = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies2 ,y_axisCol = ScaledPlotDF$PutNeuSIFT_CountAlleles, y_axisTitle = "Count derived neutral alleles")
scaledCountAllele_PutDel = plotFxn(dataFrame = ScaledPlotDF, x_axisCol = ScaledPlotDF$Subspecies2 ,y_axisCol = ScaledPlotDF$PutDelSIFT_CountAlleles, y_axisTitle = "Count derived deleterious alleles")

#Neutral
scaledPutNeu = ggarrange( scaledCountDerHom_PutNeu + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count Homozygotes"),
                          scaledCountVar_PutNeu + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count of Variants"),
                          scaledCountAllele_PutNeu + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + ggtitle("Count of Alleles"),
                    align = 'hv',
                    labels = c("A", "B", "C"),
                    hjust = -1,
                    nrow = 1)
scaledPutNeuAnnot = annotate_figure(scaledPutNeu, left = text_grob("Neutral", size=24, face="bold",rot = 90, hjust = 0.5))

#Deleterious
scaledPutDel = ggarrange( scaledCountDerHom_PutDel + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                          scaledCountVar_PutDel + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                          scaledCountAllele_PutDel + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
                    align = 'hv',
                    labels = c("D", "E", "F"),
                    hjust = -1,
                    nrow = 1)
scaledPutDelAnnot = annotate_figure(scaledPutDel, left = text_grob("Deleterious", size=24, face="bold",rot = 90, hjust = 0.5))

#plot deleterious and neutral
ggarrange(scaledPutNeuAnnot, scaledPutDelAnnot, nrow = 2)

#Compute pvalues
pairwise.wilcox.test(ScaledPlotDF$PutDelSIFT_CountDerHom, ScaledPlotDF$Subspecies2, p.adj = "bonf")$p.value
pairwise.wilcox.test(ScaledPlotDF$PutDelSIFT_CountVariants, ScaledPlotDF$Subspecies2, p.adj = "bonf")$p.value
pairwise.wilcox.test(ScaledPlotDF$PutDelSIFT_CountAlleles, ScaledPlotDF$Subspecies2, p.adj = "bonf")$p.value
  