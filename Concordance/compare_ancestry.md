##Plotting ancestry concordance
library(ggplot2)
setwd('~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Ellie’s MacBook Pro (2)/Captives-2022-Final/Concordance/')

df <- read.csv('ancestry_concordance_all.csv')

library(tidyverse)
library(RColorBrewer)

cbPalette = c("Full" = "#E69F00",  "5x" = "#56B4E9", "2x" = "#009E73", "1x" = "#F0E442", "0.5x" = "#0072B2", "0.25x"="#D55E00")#palette

pivotData = df %>% 
  pivot_longer(cols = Bengal:Indochinese, 
                     names_to = "Species", 
                     values_to = "Ancestry")
ggplot(pivotData, aes(x=Species, y=Ancestry, colour=Coverage)) + 
  geom_jitter(height = 0, width = 0.1, show.legend = FALSE) +
  scale_colour_manual(values = cbPalette) +
  labs(x="Source", y="Ancestry fraction") +
  facet_grid(Individual~.) +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 16, angle=90), 
        axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 14)) +
  theme(legend.position = "right")

Amur  <- pivotData %>%
  filter(Species == 'Amur')

Bengal  <- pivotData %>%
  filter(Species == 'Bengal')

Malayan  <- pivotData %>%
  filter(Species == 'Malayan')

Indochinese  <- pivotData %>%
  filter(Species == 'Indochinese')

Sumatran  <- pivotData %>%
  filter(Species == 'Sumatran')

SouthChina  <- pivotData %>%
  filter(Species == 'SouthChina')
  
ggplot(data = Amur, aes(x = Ancestry, fill = Coverage)) + 
  geom_density(alpha = 0.3)

ggplot(data = Bengal, aes(x = Ancestry, fill = Coverage)) + 
  geom_density(alpha = 0.3)

ggplot(data = Indochinese, aes(x = Ancestry, fill = Coverage)) + 
  geom_density(alpha = 0.3)

ggplot(data = Malayan, aes(x = Ancestry, fill = Coverage)) + 
  geom_density(alpha = 0.3)

ggplot(data = Sumatran, aes(x = Ancestry, fill = Coverage)) + 
  geom_density(alpha = 0.3)

ggplot(data = SouthChina, aes(x = Ancestry, fill = Coverage)) + 
  geom_density(alpha = 0.3)

#Pivot wider for KS test
Amurpivot = Amur %>% 
  pivot_wider(names_from = Coverage, 
               values_from = Ancestry)
Amurpivot <- Amurpivot %>%
  select(-c('Species','Individual'))

Bengalpivot = Bengal %>% 
  pivot_wider(names_from = Coverage, 
              values_from = Ancestry)
Bengalpivot <- Bengalpivot %>%
  select(-c('Species','Individual'))

Malayanpivot = Malayan %>% 
  pivot_wider(names_from = Coverage, 
              values_from = Ancestry)
Malayanpivot <- Malayanpivot %>%
  select(-c('Species','Individual'))

Indochinesepivot = Indochinese %>% 
  pivot_wider(names_from = Coverage, 
              values_from = Ancestry)
Indochinesepivot <- Indochinesepivot %>%
  select(-c('Species','Individual'))

Sumatranpivot = Sumatran %>% 
  pivot_wider(names_from = Coverage, 
              values_from = Ancestry)
Sumatranpivot <- Sumatranpivot %>%
  select(-c('Species','Individual'))

SouthChinapivot = SouthChina %>% 
  pivot_wider(names_from = Coverage, 
              values_from = Ancestry)
SouthChinapivot <- SouthChinapivot %>%
  select(-c('Species','Individual'))


# Amur KS Loop ------------------------------------------------------------

# Assume df is your data frame and base_column is the base column name
base_column <- Amurpivot$Full
# Create an empty data frame to store the results
Amurresults <- data.frame(column_name = character(), p_value = numeric(), D_value = numeric(), stringsAsFactors = FALSE)
# Loop through each column in the data frame
for (col in names(Amurpivot)) {
  if (col != "base_column") {
    # Perform the KS test and store the results in the data frame
    test_result <- ks.test(Amurpivot[[col]], base_column)
    Amurresults <- rbind(Amurresults, data.frame(column_name = col, p_value = test_result$p.value, D_value = test_result$statistic))
  }
}


# Bengal KS Loop ----------------------------------------------------------
# Assume df is your data frame and base_column is the base column name
base_column <- Bengalpivot$Full
# Create an empty data frame to store the results
Bengalresults <- data.frame(column_name = character(), p_value = numeric(), D_value = numeric(), stringsAsFactors = FALSE)
# Loop through each column in the data frame
for (col in names(Bengalpivot)) {
  if (col != "base_column") {
    # Perform the KS test and store the results in the data frame
    test_result <- ks.test(Bengalpivot[[col]], base_column)
    Bengalresults <- rbind(Bengalresults, data.frame(column_name = col, p_value = test_result$p.value, D_value = test_result$statistic))
  }
}


# IndoC KS Loop -----------------------------------------------------------
# Assume df is your data frame and base_column is the base column name
base_column <- Indochinesepivot$Full
# Create an empty data frame to store the results
Indochineseresults <- data.frame(column_name = character(), p_value = numeric(), D_value = numeric(), stringsAsFactors = FALSE)
# Loop through each column in the data frame
for (col in names(Indochinesepivot)) {
  if (col != "base_column") {
    # Perform the KS test and store the results in the data frame
    test_result <- ks.test(Indochinesepivot[[col]], base_column)
    Indochineseresults <- rbind(Indochineseresults, data.frame(column_name = col, p_value = test_result$p.value, D_value = test_result$statistic))
  }
}


# Malayan KS Loop ---------------------------------------------------------
base_column <- Malayanpivot$Full
# Create an empty data frame to store the results
Malayanresults <- data.frame(column_name = character(), p_value = numeric(), D_value = numeric(), stringsAsFactors = FALSE)
# Loop through each column in the data frame
for (col in names(Malayanpivot)) {
  if (col != "base_column") {
    # Perform the KS test and store the results in the data frame
    test_result <- ks.test(Malayanpivot[[col]], base_column)
    Malayanresults <- rbind(Malayanresults, data.frame(column_name = col, p_value = test_result$p.value, D_value = test_result$statistic))
  }
}


# Sumatran KS Loop --------------------------------------------------------
base_column <- Sumatranpivot$Full
# Create an empty data frame to store the results
Sumatranresults <- data.frame(column_name = character(), p_value = numeric(), D_value = numeric(), stringsAsFactors = FALSE)
# Loop through each column in the data frame
for (col in names(Sumatranpivot)) {
  if (col != "base_column") {
    # Perform the KS test and store the results in the data frame
    test_result <- ks.test(Sumatranpivot[[col]], base_column)
    Sumatranresults <- rbind(Sumatranresults, data.frame(column_name = col, p_value = test_result$p.value, D_value = test_result$statistic))
  }
}

# South China KS Loop -----------------------------------------------------
base_column <- SouthChinapivot$Full
# Create an empty data frame to store the results
SouthChinaresults <- data.frame(column_name = character(), p_value = numeric(), D_value = numeric(), stringsAsFactors = FALSE)
# Loop through each column in the data frame
for (col in names(SouthChinapivot)) {
  if (col != "base_column") {
    # Perform the KS test and store the results in the data frame
    test_result <- ks.test(SouthChinapivot[[col]], base_column)
    SouthChinaresults <- rbind(SouthChinaresults, data.frame(column_name = col, p_value = test_result$p.value, D_value = test_result$statistic))
  }
}

# View results ------------------------------------------------------------
print(Amurresults)
print(Bengalresults)
print(Malayanresults)
print(Indochineseresults)
print(Sumatranresults)
print(SouthChinaresults)
