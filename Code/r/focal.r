# Modules
library(readr)
library(ggplot2)
library(reshape2)
library(gghighlight)
options(scipen = 999) # prevent scientific notation for genomic position on x axis

population <- read_csv('/home/baron/Documents/PhD/ML_balancingselection_data/abi_raw_data/ngs_output/plot_data/final/region_data_2.csv',
                   col_types = cols(Position = col_character(), Allele_frequency = col_number(),
                                    Generation = col_number()))
focal <- population$Position['18521001']
#View file if required
#View(allele)

# General plot for all variants within focal region
ggplot(population, aes(x = Generation, y = Allele_frequency, color = Position)) +
  #geom_count(size = 3) +
  geom_line(linewidth = 0.5) +
  #geom_line(data = population['Position'[18521001]]) +
  #geom_smooth(method = 'lm', se = FALSE) +
  #scale_color_discrete(breaks = c(2,4,8,12,20,28,36,44,56)) +
  theme_grey() +
  #scale_color_gradient2(low = 'black', mid = 'yellow', high = 'blue')
  labs(title = 'Alternate allele frequency within focal region (18520501-18521500) - Cage 1', y = 'Alternate allele frequency', x = 'Generation') +
  scale_x_continuous(breaks = c(2,4,8,12,20,28,36,44,56)) +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))
