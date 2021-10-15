####Load libraries and set working directory
library(tidyverse)
options(scipen = 999) #use this to turn off scientific notation
setwd("~/TigerProject/FSNP_Het/")

####Load files
heterozygosity = read_delim('~/TigerProject/FSNP_Het/highcov-lowcov-nodups.ba-AN-MM-pcc-GM-MASTER.vcftoolshet.het', delim = '\t', col_names = TRUE)
names(heterozygosity) = c('Individual','O_HOM', 'E_HOM', 'N_SITES','FSNP') #rename columns

popsDF = read_csv("~/TigerProject/IndivFiles/individual_ids.csv") %>%
  mutate(combo = ifelse(Subspecies == "Unknown", paste(Subspecies, "-", Phenotype, sep = ""), Subspecies))

####Make plotting data frame 
all_nodups_full_highCov = heterozygosity %>%
  left_join(popsDF) %>%
  mutate(Heterozygosity = (N_SITES-O_HOM)/1476111759,
         Subspecies2 = factor(Subspecies2, levels = c('Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran','Unknown'))) %>%
  arrange(Subspecies2, Heterozygosity) %>%
  filter(Individual != 'T1' & Coverage == "High") #remove T1 and low coverage individuals


####Plot heterozygosity all-nodups ------------------------------------------
cbPalette = c("Generic" = "gray25", "Generic* (Amur)" = "lightsalmon2" , "Generic* (Bengal)" = "lightskyblue2", "Amur" = "#D55E00",  "Bengal" = "steelblue", "Malayan" = "#009E73", "Sumatran" = "gold3", "Indochinese" = "chocolate1", "South China" = "firebrick2", "Unknown"="darkmagenta", "Unknown-Orange" = "#CC79A7", "Unknown-SnowWhite" = "#867BCF", "Unknown-Golden"="darkseagreen3", "Unknown-White"="cornflowerblue")#palette

#order individuals
IDs = all_nodups_full_highCov$Individual 
all_nodups_full_highCov$Individual = factor(all_nodups_full_highCov$Individual, levels = IDs)

#lollipop plot
LolliPlotHet =ggplot(all_nodups_full_highCov, aes(x=Individual, y=Heterozygosity, colour = Subspecies2)) + 
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
#jpeg(file="heterozygosity-vcftools-all-nodups-high-cov.jpeg", width = 800, height = 500, res = 100)
VioPlotHet = ggplot(all_nodups_full_highCov, aes(x=Subspecies2, y=Heterozygosity)) + 
  geom_violin() +
  geom_jitter(height = 0, width = 0.1, aes(colour = Subspecies2)) +
  #geom_jitter(height = 0, width = 0.1, aes(shape = Coverage)) + #use if you leave low coverage data
  #scale_shape_manual(values=c(1,18)) + #use if you leave low coverage data
  scale_color_manual(name="Subspecies", values=cbPalette) +
  labs(x="Subspecies", y="Observed Heterozygosity") +
  theme_bw() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y = element_text(size = 18), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.position = "none")
#dev.off()


####Plot FSNP
FSNP = ggplot(all_nodups_full_highCov, aes(x=Subspecies2, y=FSNP)) +
  geom_violin(aes(fill=Subspecies2)) + 
  geom_jitter(height = 0, width = 0.1) +
  scale_fill_manual(name = "Subspecies", values = cbPalette) + 
  labs(x="Subspecies", y=expression(F[SNP])) + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))
