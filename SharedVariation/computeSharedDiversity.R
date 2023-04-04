####Load libraries and set working directory
library(tidyverse)
setwd("~/Documents/Tigers/SharedVariation/")

####Functions


####Read in all the wild tigers
df = read_delim("combineDataReps.counts") %>%
  select(!c("CHROM", "POS")) %>%
  mutate_all(function(x) gsub("C:|T:|A:|G:","",x)) %>%
  mutate_if(is.character, as.numeric)

####Iterate through replicates to get the fraction of segregating sites that are fixed in wild tigers

reps = colnames(df)[3:12] #columns with replicates
segSitesOnly = data.frame()

for (i in reps) {
  prepData = df %>%
    select(i, Amur, Bengal, Indoc, Malayan, Sumatran) %>% #only variant sites in captives
    mutate(atLeastOneSNPWild = ifelse(Amur==0&Bengal==0&Indoc==0&Malayan==0&Sumatran==0, 0, 1))
  
  #rename columns
  colnames(prepData) = c("Generic", "Amur", "Bengal", "Indochinese", "Malayan", "Sumatran", "atLeastOneSNPWild")
  
  removeFixedSitesGen = prepData %>%
    filter(Generic != 0)
    
    countsDF = colSums(removeFixedSitesGen != 0) %>% 
    as.data.frame() %>% 
    rownames_to_column(var ="Subspecies") %>%
    rename("segSites"=".") %>%
    pivot_wider(names_from = Subspecies,
                values_from = segSites)  %>% 
    mutate_at(vars(-matches("Generic")), ~ . -Generic) %>% #get the number of fixed sites
    mutate_if(is.numeric, abs) %>% 
    mutate_at(vars(-matches("Generic")), ~ . /Generic) %>% #proportion 
    rename("fixedSites"="atLeastOneSNPWild")
  
  segSitesOnly = rbind(segSitesOnly, countsDF)
}
row.names(segSitesOnly) = seq(1:10)

summary_stats = segSitesOnly %>%
  select(-c("Generic")) %>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>%
  group_by(name) %>%
  summarise(mean = mean(value, na.rm=TRUE),
            sd = sd(value, na.rm=TRUE)) %>% 
  mutate(name = gsub("fixedSites", "Fixed in wild", name),
         Subspecies2 = factor(name, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran', 'Fixed in wild')))

#plot 
cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Indochinese" = "gold4", "South China" = "plum", "Sumatran" = "cornflowerblue", "Fixed in wild"="black")#palette

ggplot(summary_stats, aes(x=Subspecies2, y = mean, colour= Subspecies2)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, size=2) + 
  geom_point(color="black", size=3) +
  #coord_flip() +  
  scale_colour_manual(name = "Subspecies", values = cbPalette) +
  labs(x = "", y="Proportion of fixed sites\nthat are segregating in captives") +
  theme_bw() + 
  theme(axis.text.x = element_text( hjust= 0.5, vjust=1, size=20), 
        axis.text.y = element_text(size = 20), 
        plot.title=element_text(size =24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=24), 
        legend.position = "none") 
