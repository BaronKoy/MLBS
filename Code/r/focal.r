# Modules
library(readr)
library(ggplot2)
library(reshape2)
options(scipen = 999) # prevent scientific notation for genomic position on x axis

population <- read_csv('/home/baron/Documents/PhD/ML_balancingselection_data/reuter_data/abi_pileups/focal_outputs/plot_data/region_data_1.csv',
                   col_types = cols(Position = col_character(), Allele_frequency = col_number(),
                                    Generation = col_number()))
#View file if required
#View(allele)

cage_plot <- ggplot(population, aes(x = Generation, y = Allele_frequency, color = Position)) +
  #geom_count(size = 3) +
  geom_line(size = 0.5) +
  #geom_line(data = population['Position'[18521001]]) +
  #geom_smooth(method = 'lm', se = FALSE) +
  #scale_color_discrete(breaks = c(2,4,8,12,20,28,36,44,56)) +
  theme_grey() +
  labs(title = 'Alternate allele frequency within focal region (18520501-18521500) - Cage 1', y = 'Alternate allele frequency', x = 'Generation') +
  scale_x_continuous(breaks = c(2,4,8,12,20,28,36,44,56)) +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))
  
cage_plot