# Modules
library(readr)
library(ggplot2)
library(reshape2)
options(scipen = 999) # prevent scientific notation for genomic position on x axis

freq <- read_csv('/home/baron/Documents/PhD/ML_balancingselection/reuter_data/pileup/depth/output/pop_10/final_pop10.pileup', #change for specific alignment file
                 col_types = cols(pos = col_number(), depth = col_number(),
                                  Generation = col_character()))
#assign column types

#View file if required
#View(freq)

line_plot <- ggplot(freq, aes(x = pos, y = depth, color = Generation)) +
  geom_line(size = 0.75) +
  theme_grey() +
  labs(title = 'Distribution depth across focal region - pop_10, gen2-56', y = 'Read_depth(n=96)', x = 'Genomic_position_(Chr_3R)') +
  #scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,96)) + 
  #scale_y_continuous(breaks = c(0,10000,20000,3000,40000, 50000)) +
  #scale_fill_continuous(breaks = c(0,2,4,8,12,20,28,36,44,56)) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30), 
        legend.key.size = unit(1.7, 'cm')) + theme(plot.title = element_text(size = 30)) +
  scale_color_discrete(breaks = c(0,2,4,8,12,20,28,36,44,56))

line_plot