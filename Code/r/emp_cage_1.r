# Modules
library(readr)
library(ggplot2)
library(reshape2)
options(scipen = 999) # prevent scientific notation for genomic position on x axis

population <- read_csv('/home/baron/Documents/PhD/ML_balancingselection_data/reuter_data/abi_pileups/focal_outputs/columns/cage_1.csv',
                   col_types = cols(Position = col_number(), Allele_frequency = col_number(),
                                    Generation = col_character()))
#View file if required
#View(allele)

cage_plot <- ggplot(population, aes(x = Position, y = Allele_frequency, color = Generation)) +
  geom_point(size = 3) + 
  geom_line(size = 1) +
  #geom_smooth(method = 'lm', se = FALSE) +
  labs(title = 'Allele frequency per generation (Cage 1)', y = 'Alternate allele frequency', x = 'Genomic position (chromosome 3R)') +
  scale_color_discrete(breaks = c(2,4,8,12,20,28,36,44,56)) + theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15)) +
  theme_grey()
cage_plot