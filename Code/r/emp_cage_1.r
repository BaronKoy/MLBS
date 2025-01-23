# Modules
library(readr)
library(ggplot2)
library(reshape2)
options(scipen = 999) # prevent scientific notation for genomic position on x axis

population <- read_csv('/home/baron/Documents/PhD/ML_balancingselection/reuter_data/abi_pileups/outputs/plot_data/cage_1.csv',
                   col_types = cols(Position = col_number(), Allele_Frequency = col_number(),
                                    Generation = col_character()))
#View file if required
#View(allele)

cage_plot <- ggplot(population, aes(x = Position, y = Allele_Frequency, color = Generation))+
  geom_line() +
  theme_grey() +
  labs(title = 'Allele frequency per generation (Cage 1)', y = 'Alternate allele frequency', x = 'Genomic position (chromosome 3R)') +
  scale_color_discrete(breaks = c(2,4,8,12,20,28,36,44,56)) + theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))
  
cage_plot