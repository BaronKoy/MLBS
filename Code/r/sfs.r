# Modules
library(readr)
library(ggplot2)
library(reshape2)
options(scipen = 999) # prevent scientific notation for genomic position on x axis

freq <- read_csv('/home/baron/Documents/rotation_2/QM_rotation/scripts/outputs/vcfs/cage_1/neutral/all.csv',
                   col_types = cols(Frequency = col_number(), Generation = col_number()))
#View file if required
#View(allele)

cage_plot <- ggplot(freq, aes(x = Frequency, fill = Generation, group = Generation)) + geom_histogram(bins = 96) +
  theme_grey() +
  labs(title = 'Site frequency spectrum for neutrality', y = 'Proportion of Alt allele', x = 'Alt allele frequency class', color = 'Generations') + 
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,96)) + 
  scale_y_continuous(breaks = c(0,10000,20000,3000,40000, 50000)) +
  scale_fill_continuous(breaks = c(0,2,4,8,12,20,28,36,44,56)) +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15), legend.key.size = unit(1.7, 'cm'))
  
cage_plot