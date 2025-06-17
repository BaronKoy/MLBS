# Modules
library(readr)
library(ggplot2)
library(reshape2)
options(scipen = 999) # prevent scientific notation for genomic position on x axis

freq <- read_csv('/home/baron/Documents/PhD/Data/abi_raw_data/BAMs/cage_1/final.txt', #change "cage_1" string for specific alignment file
                 col_types = cols(pos = col_number(), mapping_quality = col_number(),
                                  Generation = col_character()))
#assign column types

#View file if required
#View(freq)

line_plot <- ggplot(freq, aes(x = pos, y = mapping_quality, color = Generation)) +
  geom_line(linewidth = 0.75) +
  theme_grey() +
  labs(title = 'Mapping quality across focal region, gen2-56 (Cage 1)', y = 'Mapping_quality', x = 'Genomic_position_(Chr_3R)') + # Change cage number in title as required
  #scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,96)) + 
  #scale_y_continuous(breaks = c(0,10000,20000,3000,40000, 50000)) +
  #scale_fill_continuous(breaks = c(0,2,4,8,12,20,28,36,44,56)) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30), 
        legend.key.size = unit(1.7, 'cm')) + theme(plot.title = element_text(size = 30)) +
  scale_color_discrete(breaks = c(2,4,8,12,20,28,36,44,56))

line_plot
