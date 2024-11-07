# Modules
library(readr)
library(ggplot2)
library(reshape2)
options(scipen = 999) # prevent scientific notation for genomic position on x axis

freq <- read_csv('/home/baron/Documents/PhD/ML_balancingselection/reuter_data/pileup/depth/output/S010817.pileup', #change for specific alignment file
                 col_types = cols(pos = col_number(), depth = col_number()))
#assign column types

#View file if required
#View(freq)

line_plot <- ggplot(freq, aes(x = pos, y = depth)) + geom_line(color = 'steelblue', size = 0.75) +
  theme_grey() +
  labs(title = 'Distribution depth across focal region - pop_1, gen2', y = 'Read_depth(/96)', x = 'Genomic_position_(Chr_3R)') +
  #scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,96)) + 
  #scale_y_continuous(breaks = c(0,10000,20000,3000,40000, 50000)) +
  #scale_fill_continuous(breaks = c(0,2,4,8,12,20,28,36,44,56)) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30), 
        legend.key.size = unit(1.7, 'cm')) + theme(plot.title = element_text(size = 30))

line_plot