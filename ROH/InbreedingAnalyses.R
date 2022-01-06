#Load libraries and set working directory
library(tidyverse)
library(data.table)
library(ggpubr)
setwd("~/TigerProject/ROH/")

#Load roh, FSNP, deleterious vars, and meta data
popsDF = read_csv("~/TigerProject/IndivFiles/individual_ids.csv")

heterozygosity = read_delim('~/TigerProject/FSNP_Het/highcov-lowcov-nodups.ba-AN-MM-pcc-GM-MASTER.vcftoolshet.het', delim = '\t', col_names = TRUE)
names(heterozygosity) = c('INDV','O_HOM', 'E_HOM', 'N_SITES','FSNP') #rename columns

####Annotations 
annotations = read_delim("~/TigerProject/AnnotSites/GTAnnotationCountResults_Sept2021_Tigers.txt", delim = "\t") %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  mutate(CallableSites = LineCount - Missing)


PropAnnotations = annotations[,c(1:12, 21:26)] %>%
  mutate_at(vars(LOF_CountAlleles: PutNeuSIFT_CountVariants), list(~./annotations$CallableSites)) #Just run replace on plotting fxns with regular df to make everything proportional

#merge roh and heterozygosity
roh = read_delim("~/TigerProject/ROH/TrueROH_propCoveredwithin1SDMean_gr100kb_allChroms_highCov_runSpeciesSep_garlic.txt", delim = "\t") %>%
  left_join(popsDF, by = c("INDV" = "Individual")) %>%
  left_join(heterozygosity) %>%
  mutate(Heterozygosity = (N_SITES-O_HOM)/1476111759,
         Subspecies2 = factor(Subspecies2, levels = c('Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran','Generic')))

####
rohLengthsClass = roh %>%
  group_by(INDV, TYPE) %>%
  summarise_at(c("AUTO_LEN"), sum) %>%
  ungroup() %>%
  left_join(popsDF,by = c("INDV" = "Individual"))

FROH = rohLengthsClass %>%
  filter(TYPE == "C") %>%
  mutate(Froh = AUTO_LEN/2174711735,
         Fsnp = roh$FSNP[match(INDV,roh$INDV)],
         DelHom = PropAnnotations$PutDelSIFT_CountDerHom[match(INDV,PropAnnotations$ID)],
         DelAllele = PropAnnotations$PutDelSIFT_CountAlleles[match(INDV,PropAnnotations$ID)],
         DelVar = PropAnnotations$PutDelSIFT_CountVariants[match(INDV,PropAnnotations$ID)]) 

compareMeasures = FROH %>%
  select(Subspecies2, Froh, Fsnp) %>%
  pivot_longer(!Subspecies2, names_to = "Measure", values_to = "value") %>%
  mutate(Subspecies2 = factor(Subspecies2, levels = c('Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran','Generic')))


ggplot(compareMeasures, aes(x=Measure, y=value)) +
  geom_violin() +
  facet_wrap(~Subspecies2, nrow = 2) +
  geom_violin(aes(fill=Subspecies2)) + 
  geom_jitter(height = 0, width = 0.1) +
  scale_fill_manual(name = "Subspecies", values = cbPalette) + 
  labs(x="Inbreeding Measure", y="Coefficient") + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14))
