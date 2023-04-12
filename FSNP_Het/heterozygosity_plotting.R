#heterozygosity
setwd('~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Ellie’s MacBook Pro (2)/Captives-2022-Final/Heterozygosity/')
library(ggplot2)
library(dplyr)

#plot all
heterozygosity <- read.csv('highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.lowcov86.filter.nodups.het', header = TRUE, sep='\t')
metadata <- read.table('~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Ellie’s MacBook Pro (2)/captives-new-2021/Metadata-files/plotting_metadata.csv', header=TRUE, sep = ',',quote="")

missing <- read.csv('highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.lowcov86.filter.nodups.imiss', header = TRUE, sep='\t')

missing <- missing %>%
  select(-c(N_DATA,N_GENOTYPES_FILTERED)) %>%
  dplyr::rename(Individual = INDV)
heterozygosity <- heterozygosity %>%
  dplyr::rename(Individual = INDV) %>%
  mutate(heterozygosity=((N_SITES-O.HOM.)/2174711735))

het_metadata <- merge(heterozygosity, metadata, by="Individual")
het_meta_missing <- merge(het_metadata, missing, by="Individual")

cbPalette = c("Amur" = "#0072B2",  "Bengal" = "#882255", "Malayan" = "#009E73", "Sumatran" = "cornflowerblue", "Indochinese" = "gold4", "Generic"="gray25", "SouthChina" = "plum")#palette

#plot with imputed identified
ggplot(het_meta_missing, aes(x=Subspecies_GroupID_2, y=heterozygosity, fill=Subspecies_GroupID_2)) + 
  geom_violin(position=position_dodge(1)) + 
  scale_fill_manual(values = cbPalette, name = "") +
  geom_point(aes(shape=Coverage_group), position=position_jitter(0.2), size = 0.6) +
  scale_shape_manual(values=c(1, 2), name = "") +
  ylab('Observed Heterozygosity') + 
  xlab('') + 
  theme_bw() +
  theme(axis.title=element_text(size = 14), axis.text.x=element_text(size = 12)) +
  labs(color=NULL)


#correlate to missingness
het_meta_missing_noimpute <- het_meta_missing %>%
  filter(Coverage_group == 'Unimputed')

ggplot(het_meta_missing_noimpute, aes(x = F_MISS, y = heterozygosity)) +
  geom_point(color = "#00AFBB", size = 2, shape = 19) +
  xlab('Missing') +
  ylab('het') + 
  theme_bw()
  
#heterozygosity with no samples >F_MISS 0.1
het_meta_missing_noimpute_FMiss <- het_meta_missing_noimpute %>%
  filter(F_MISS <= 0.20)

het_meta_missing_noimpute_FMiss %>%
  group_by(Subspecies_GroupID_2) %>%
  summarize(min = min(heterozygosity),
            q1 = quantile(heterozygosity, 0.25),
            median = median(heterozygosity),
            mean = mean(heterozygosity),
            q3 = quantile(heterozygosity, 0.75),
            max = max(heterozygosity))

tapply(het_meta_missing_noimpute_FMiss$heterozygosity,het_meta_missing_noimpute_FMiss$Subspecies_GroupID_2,FUN = var)

ggplot(het_meta_missing_noimpute_FMiss, aes(x=Subspecies_GroupID_2, y=heterozygosity, fill=Subspecies_GroupID_2)) + 
  geom_violin(position=position_dodge(1)) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  ylab('Observed Heterozygosity') + 
  xlab('') + 
  theme_bw() 

