#Load libraries and data
library(tidyverse)
library(data.table)
library(ggpubr)

#Load ibdseq roh outputs and meta data
setwd("~/TigerProject/ROH/")

popsDF = read_csv("~/TigerProject/IndivFiles/individual_ids.csv")

roh = read_delim("TrueROH_propCoveredwithin1SDMean_gr100kb_allChroms_highCov_runSpeciesSep_garlic.txt", delim = "\t") %>%
  left_join(popsDF,by = c("INDV" = "Individual"))

rohLengthsClass = roh %>%
  group_by(INDV, TYPE) %>%
  summarise_at(c("AUTO_LEN"), sum) %>%
  ungroup() %>%
  left_join(popsDF,by = c("INDV" = "Individual"))

FROH = rohLengthsClass %>%
  filter(TYPE == "C") %>%
  mutate(Froh = AUTO_LEN/2174711735)
summaryTable = roh %>%
  arrange(TYPE, Subspecies) %>%
  group_by(TYPE, Subspecies) %>%
  summarise(mean = mean(AUTO_LEN)/10^3, 
            min = min(AUTO_LEN)/10^3, 
            max = max(AUTO_LEN)/10^3) %>%
  mutate_if(is.numeric, round, digits=3) %>%
  ungroup() %>%
  rename("ROH class" = TYPE, "subspecies" = "Subspecies")

##Plot data
cbPalette = c("Generic" = "gray25", "Generic* (Amur)" = "lightsalmon2" , "Generic* (Bengal)" = "lightskyblue2", "Amur" = "#D55E00",  "Bengal" = "steelblue", "Malayan" = "#009E73", "Sumatran" = "gold3", "Indochinese" = "chocolate1", "South China" = "firebrick2", "Unknown"="darkmagenta", "Unknown-Orange" = "#CC79A7", "Unknown-SnowWhite" = "#867BCF", "Unknown-Golden"="darkseagreen3", "Unknown-White"="cornflowerblue") #palette


ggtexttable(summaryTable, rows = NULL, theme = ttheme("mBlackWhite"))

ggplot(rohLengthsClass, aes(x=Subspecies, y=AUTO_LEN/10^6, fill=Subspecies)) + 
  geom_boxplot() +
  facet_wrap(~TYPE, scale="free") +
  scale_fill_manual(name = "Sub-species", values = cbPalette) +
  labs(x="Sub-species", y="Length of genome in ROH (Mb)") +
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14))  

ggplot(FROH, aes(x=Subspecies, y=Froh, colour=Subspecies)) +
  geom_boxplot(size=1) + 
  geom_point(size=0.5) + 
  scale_colour_manual(name = "Sub-species", values = cbPalette) + 
  labs(x = "Sub-species", y=expression(F[ROH])) + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))

ggplot(FROH, aes(x=Subspecies, y=Froh, colour=Subspecies)) +
  geom_boxplot(size=1) + 
  geom_point(size=0.5) + 
  scale_colour_manual(name = "Sub-species", values = cbPalette) + 
  labs(x = "Sub-species", y=expression(F[ROH])) + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))