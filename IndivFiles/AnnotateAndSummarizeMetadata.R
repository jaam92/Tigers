#Libraries and functions
load(tidyverse)
`%ni%` <- Negate(`%in%`)

unrelateds = read_delim("~/Documents/Tigers/Relatedness/FindDups/SNPRelate/unrelateds_pcair/N10-N6-N3_unrelateds.txt") %>% 
  mutate(Subspecies2 = NULL)
colnames(unrelateds) = c("Individual", "SFS_N3", "SFS_N6", "SFS_N10")

removed = c('T18', 'T5', 'T10', 'SRR5591010', 'SRR5612311', 'SRR5612312')

#Make annotated file
popsDF = read_csv("~/Documents/Tigers/IndivFiles/individual_ids.csv") %>%
  mutate(Subspecies2 = factor(Subspecies2, levels = c('Generic', 'Amur', 'Bengal', 'Indochinese', 'Malayan', 'South China', 'Sumatran'))) %>%
  select(Individual, Subspecies2, Coverage, Phenotype, Sex) %>%
  mutate(PCA=ifelse(Individual%ni%removed, 1, 0),
         VCFTools_Relatedness = 1,
         TRUFFLE_IBDhalf = 1,
         Pedigree_Kinship = 0,
         PC_air = ifelse(Individual%ni%removed, 1, 0),
         Unrelated_group = ifelse(Individual%in%unrelateds$Individual, 1, 0),
         Sex_detection = ifelse(Individual%ni%removed, 1, 0),
         Mitochrondrial_haplotypes = ifelse(Individual%ni%removed, 1, 0),
         Heterozygosity = ifelse(Individual%ni%removed & Coverage == "High", 1, 0),
         Runs_of_homozygosity = ifelse(Individual%ni%removed & Coverage == "High", 1, 0),
         Identity_by_descent =  ifelse(Individual%ni%removed & Coverage == "High", 1, 0),
         Variant_annotation = ifelse(Individual%ni%removed & Coverage == "High", 1, 0)) %>%
  rename(Subspecies = Subspecies2) %>%
  left_join(unrelateds) %>%
  mutate(Reason_for_exclusion = ifelse(Individual%in%removed, "duplicate samples", ""),
         Reason_for_exclusion = ifelse(Individual%ni%removed & Coverage == "Low", "low coverage", Reason_for_exclusion)) 

popsDF[is.na(popsDF)] <- 0 #replace na with 0

#write.csv(popsDF, file = "~/Documents/Tigers/IndivFiles/individual_ids_annotated.csv") write it out

#All the counts
start = popsDF %>%
  group_by(Subspecies2) %>%
  count()


filtered = popsDF %>%
  filter(!Individual %in% removed) 

subspecies = filtered %>%
  group_by(Subspecies2) %>%
  count(name = "N")
colnames(subspecies) = c("Subspecies", "Sample size")

coverage = filtered %>%
  group_by(Subspecies2, Coverage) %>%
  count(name = "N")
colnames(coverage) = c("Subspecies", "Coverage", "Sample size")

phenotypeCov = filtered %>%
  mutate(Subpheno = paste(Subspecies2, Phenotype, sep = "-")) %>%
  group_by(Subpheno, Coverage) %>%
  count(name = "N")
colnames(phenotypeCov) = c("Subspecies (with phenotype)", "Coverage", "Sample size")  

plotSubspecies = ggtexttable(subspecies, rows = NULL, theme = ttheme("mBlackWhite"))
plotCov = ggtexttable(coverage, rows = NULL, theme = ttheme("mBlackWhite"))
plotPheno = ggtexttable(phenotypeCov, rows = NULL, theme = ttheme("mBlackWhite"))

ggarrange(plotSubspecies, plotCov, plotPheno, ncol = 3, labels = c("A","B", "C"), font.label = list(size = 18), align = 'hv')
