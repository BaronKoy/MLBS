# Modules
library(readr)
library(ggplot2)
library(reshape2)
options(scipen = 999) # prevent scientific notation for genomic position on x axis

freq <- read_csv('/home/baron/Documents/PhD/ML_balancingselection/reuter_data/abi_pileups/outputs/plot_data/S010817_final.txt',
                 col_types = cols(position = col_number(), saf_MLE = col_number(),
                                  dataset = col_character()))
#assign column types

#View file if required
#View(freq)

line_plot <- ggplot(freq, aes(x = position, y = saf_MLE, color = dataset)) +
  geom_point(size = 0.75) +
  theme_grey() +
  labs(title = 'Minor allele frequency likelihoods - pop_1, gen2', y = 'saf_MLE', x = 'Genomic_position_(Chr_3R)') +
  #scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,96)) + 
  #scale_y_continuous(breaks = c(0,10000,20000,3000,40000, 50000)) +
  #scale_fill_continuous(breaks = c(0,2,4,8,12,20,28,36,44,56)) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30), 
        legend.key.size = unit(1.7, 'cm')) + theme(plot.title = element_text(size = 30))

line_plot