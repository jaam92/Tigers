#Load libraries and set working directory
library(tidyverse)
library(data.table)
library(ggpubr)
setwd("~/Documents/Tigers/ROH/")

#Load roh and FSNP and meta data
popsDF = read_csv("~/Documents/Tigers/IndivFiles/individual_ids.csv")

roh = read_delim("~/Documents/Tigers/ROH/TrueROH_propCoveredwithin1SDMean_gr100kb_allChroms_highCov_runSpeciesSep_garlic.txt", delim = "\t") %>%
  left_join(popsDF, by = c("INDV" = "Individual")) %>%
  mutate(Subspecies2 = factor(Subspecies2, levels = c('Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran','Generic')))

rohLengthsClass = roh %>%
  group_by(INDV, TYPE) %>%
  summarise_at(c("AUTO_LEN"), sum) %>%
  ungroup() %>%
  left_join(popsDF,by = c("INDV" = "Individual"))

FROH = rohLengthsClass %>%
  filter(TYPE == "C") %>%
  mutate(Froh = AUTO_LEN/2174711735)

summaryTable = roh %>%
  arrange(TYPE, Subspecies2) %>%
  group_by(TYPE, Subspecies2) %>%
  summarise(mean = mean(AUTO_LEN)/10^3, 
            min = min(AUTO_LEN)/10^3, 
            max = max(AUTO_LEN)/10^3) %>%
  mutate_if(is.numeric, round, digits=3) %>%
  ungroup() %>%
  rename("ROH class" = TYPE, "subspecies" = "Subspecies2")

ggtexttable(summaryTable, rows = NULL, theme = ttheme("mBlackWhite"))

##Plot data
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Sumatran" = "cornflowerblue", "Indochinese" = "gold4", "South China" = "plum", "Generic"="gray25")#palette
cbPalette_expanded = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Sumatran" = "cornflowerblue", "Indochinese" = "gold4", "South China" = "plum", "Generic"="gray25", "Generic-Orange" = "#CC79A7", "Generic-SnowWhite" = "#867BCF", "Generic-Golden"="darkseagreen3", "Generic-White"="cornflowerblue")#palette


plotROHs = ggplot(rohLengthsClass, aes(x=Subspecies2, y=AUTO_LEN/10^6, fill=Subspecies2)) + 
  geom_boxplot() +
  facet_wrap(~TYPE, scale="free") +
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


