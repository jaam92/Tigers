####Load libraries and set working directory
library(tidyverse)
options(scipen = 999) #use this to turn off scientific notation
setwd("~/Documents/Tigers/FSNP_Het")

####Load files
heterozygosity = read_delim('~/Documents/Tigers/FSNP_Het/allSamp.het', delim = '\t', col_names = TRUE)
colnames(heterozygosity) = c('Individual','O_HOM', 'E_HOM', 'N_SITES','FSNP') #rename columns

popsDF = read_csv("~/Documents/Tigers/IndivFiles/individual_ids.csv") %>%
  mutate(combo = ifelse(Subspecies2 == "Generic", paste(Subspecies2, "-", Phenotype, sep = ""), Subspecies2))

####Make plotting data frame 
all_nodups_full_highCov = heterozygosity %>%
  left_join(popsDF) %>%
  filter(!Individual %in% indivsToRemove) %>%
  mutate(Heterozygosity = (N_SITES-O_HOM)/1476111759,
         Subspecies2 = factor(Subspecies2, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran')),
         BoundFSNP = ifelse(FSNP < 0, 0, FSNP),
         Subspecies = NULL) %>%
  arrange(Subspecies2, Heterozygosity) 


####Plot heterozygosity all-nodups ------------------------------------------
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Indochinese" = "gold4", "Sumatran" = "cornflowerblue", "South China" = "plum", "Generic"="gray25")#palette
cbPalette_expanded = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Sumatran" = "cornflowerblue", "Indochinese" = "gold4", "South China" = "plum", "Generic"="gray25", "Generic-Orange" = "#CC79A7", "Generic-SnowWhite" = "#867BCF", "Generic-Golden"="darkseagreen3", "Generic-White"="cornflowerblue")#palette

#order individuals
IDs = all_nodups_full_highCov$Individual 
all_nodups_full_highCov$Individual = factor(all_nodups_full_highCov$Individual, levels = IDs)

#lollipop plot
LolliPlotHet = ggplot(all_nodups_full_highCov, aes(x=Individual, y=Heterozygosity, colour = Subspecies2)) + 
  geom_segment(aes(x=Individual, xend=Individual, y=0, yend=Heterozygosity), size = 1, colour = "grey")+
  geom_point(size = 2) +
  scale_color_manual(name="Subspecies", values=cbPalette) +
  labs(x="Individual", y="Observed Heterozygosity") +
  theme_bw() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y = element_text(size = 18), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18),
        panel.border = element_blank(),
        panel.grid = element_blank())

#violin plots
#jpeg(file="heterozygosity-vcftools-calledPerSubspecies-nodups-high-cov.jpeg", width = 1100, height = 550, res = 100)
VioPlotHet = ggplot(all_nodups_full_highCov, aes(x=Subspecies2, y=Heterozygosity)) +
  geom_violin(aes(fill=Subspecies2)) + 
  geom_jitter(height = 0, width = 0.1) +
  scale_fill_manual(name = "Subspecies", values = cbPalette) + 
  labs(x="Subspecies",  y="Observed Heterozygosity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))
#dev.off()


####Plot FSNP
#jpeg(file="FSNP_vioplot.jpeg", width = 1040, height = 700, res = 100)
plotFSNP = ggplot(all_nodups_full_highCov, aes(x=Subspecies2, y=FSNP)) +
  geom_violin(aes(fill=Subspecies2)) + 
  geom_jitter(height = 0, width = 0.1) +
  scale_y_continuous(breaks=seq(0,1,0.1)) +
  scale_fill_manual(name = "Subspecies", values = cbPalette) + 
  labs(x="Subspecies", y=expression(F[SNP])) + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14)) 

#dev.off()
