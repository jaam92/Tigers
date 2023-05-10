#Load libraries and set working directory
library(tidyverse)
library(data.table)
library(ggpubr)
setwd("~/Documents/Tigers/ROH/")

#Load roh and FSNP and meta data
popsDF = read_delim("~/Documents/Tigers/IndivFiles/TableX-SampleDetails.txt")

roh = read_delim("~/Documents/Tigers/ROH/TrueROH_propCoveredwithin1SDMean_gr100kb_allChroms_highCov_runSpeciesSep_garlic.txt", delim = "\t") %>%
  left_join(popsDF, by = c("INDV" = "Sample")) %>%
  mutate(Subspecies2 = factor(Subspecies_GroupID_Corrected, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran')),
         TYPE2 = paste("Type",TYPE))

rohLengthsClass = roh %>%
  group_by(INDV, TYPE) %>%
  summarise_at(c("AUTO_LEN"), sum) %>%
  ungroup() %>%
  left_join(popsDF,by = c("INDV" = "Sample")) %>%
  mutate(Subspecies2 = factor(Subspecies_GroupID_Corrected, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran')),
         TYPE2 = paste("Type",TYPE))

FROH = rohLengthsClass %>%
  filter(TYPE == "C") %>%
  mutate(Froh = AUTO_LEN/2174711735) 

summaryTable = roh %>%
  mutate(TYPE2 = paste("Type",TYPE)) %>%
  arrange(TYPE2, Subspecies2) %>%
  group_by(TYPE2, Subspecies2) %>%
  summarise(mean = mean(AUTO_LEN)/10^3, 
            min = min(AUTO_LEN)/10^3, 
            max = max(AUTO_LEN)/10^3) %>%
  mutate_if(is.numeric, round, digits=3) %>%
  ungroup() %>%
  rename("ROH class" = TYPE2, "Subspecies" = "Subspecies2") 

ggtexttable(summaryTable, rows = NULL, theme = ttheme("mBlackWhite"))

##Plot data
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Indochinese" = "gold4", "South China" = "plum", "Sumatran" = "cornflowerblue", "Generic"="gray25")#palette
cbPalette_expanded = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Indochinese" = "gold4", "South China" = "plum", "Sumatran" = "cornflowerblue", "Generic"="gray25", "Generic-Orange" = "#CC79A7", "Generic-SnowWhite" = "#867BCF", "Generic-Golden"="darkseagreen3", "Generic-White"="cornflowerblue")#palette


plotROHs = ggplot(rohLengthsClass, aes(x=Subspecies2, y=AUTO_LEN/10^6, fill=Subspecies2)) + 
  geom_violin() +
  geom_jitter(height = 0, width = 0.1, show.legend = "FALSE") +
  facet_wrap(~TYPE2) +
  scale_fill_manual(name = "Subspecies", values = cbPalette) +
  labs(x="Subspecies", y="Length of genome in ROH (Mb)") +
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14))  
print(plotROHs)

plotFROH = ggplot(FROH, aes(x=Subspecies2, y=Froh)) +
  geom_violin(aes(fill=Subspecies2)) + 
  geom_jitter(height = 0, width = 0.1) +
  scale_fill_manual(name = "Subspecies", values = cbPalette) + 
  labs(x="Subspecies", y=expression(F[ROH])) + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14 )) 

print(plotFROH)

###plot inbreeding measures
all_nodups_full_highCov = read_delim("~/Documents/Tigers/FSNP_Het/highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.lowcov86.filter.nodups.het", delim = "\t")
y = all_nodups_full_highCov %>%
  left_join(popsDF,by = c("INDV" = "Sample")) %>%
  mutate(Subspecies = factor(Subspecies_GroupID_Corrected, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran'))) %>%
  select(INDV, "F", Subspecies) %>%
  mutate(Type = "F[SNP]") %>%
  filter(INDV %in% FROH$INDV)

x = FROH %>%
  select(INDV, Froh, Subspecies2) %>%
  mutate(Type = "F[ROH]") 

inbreedingCoeff = x %>% 
  left_join(y, by = c("INDV")) %>%
  na.omit() 

inbreeding = ggplot(inbreedingCoeff, aes(x=Subspecies2, y=FROH, fill = Subspecies2)) +
  geom_violin(scale = "width", alpha=0.8) + 
  geom_point(aes(colour = F)) + 
  scale_x_discrete(limits=c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran')) +
  scale_colour_gradient(low = "black", high = "darkorange", name=expression(F[SNP])) +
  scale_fill_manual(values=cbPalette, guide = "none") + 
  labs(x= "Subspecies", y = expression(F[ROH])) + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14 )) 


figure2 = ggarrange(plotROHs, 
                 inbreeding + xlab(NULL), 
                 nrow = 2)



####Plot FSNP versus FROH linear regression color by number of deleterious homs
source(file = "~/Documents/Tigers/AnnotSites/PlotAnnotations.R")
colnames(x) = c("INDV", "FROH","Subspecies2" ,"TYPE")
colnames(y) = c("INDV", "FSNP","Subspecies2" ,"TYPE")

x %>% 
  filter(Subspecies2 != "South China") %>%
  left_join(y, by = c("INDV")) %>%
  mutate(PutDelSIFT_CountDerHom = ScaledPlotDF$PutDelSIFT_CountDerHom[match(INDV, ScaledPlotDF$ID)]) %>%
  na.omit() %>%
  ggplot(aes(x=FSNP, y=FROH)) + 
  geom_point(aes(colour = cut(PutDelSIFT_CountDerHom, c(0, 40, 60, 80, 100))), size = 3) +
  geom_smooth(method = "lm") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), color = "red", geom = "label") + #this is not the adjusted r2
  facet_wrap(~Subspecies2.x) + 
  labs(x=expression(F[SNP]), y=expression(F[ROH])) + 
  scale_color_manual(name = "Count putatively deleterious\n derived homozygotes", 
                     values = c("(0,40]" = "black",
                                "(40,60]" = "yellow", 
                                "(60,80]" = "orange", 
                                "(80,100]" = "red"),
                     labels = c("<=40",
                                "40 < variants <= 60", 
                                "60 < variants <= 80", 
                                "80 < variants <= 100"))  +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16))
