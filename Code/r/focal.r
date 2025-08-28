# Modules
library(readr)
library(ggplot2)
library(reshape2)
library(gghighlight)
library(magrittr)
library(dplyr) 

options(scipen = 999) # prevent scientific notation for genomic position on x axis

population <- read_csv('/home/baron/Documents/PhD/Data/plot_data/cage_1/cage_1.csv', # Change "Cage_1.csv" string as required
                   col_types = cols(Position = col_character(), Allele_frequency = col_number(),
                                    Generation = col_number()))
#population$Positioncolour <- rep('black', nrow(population)) # If filtering is required
#poscolour = as.character(population$Positioncolour)

#View file if required
#View(population)

# General plot for all variants within focal region
ggplot(population, aes(x = Generation, y = Allele_frequency, color = Position)) +
  # Colour data points that appear within the focal region. Note this will depend on what alleles pop up in what cage populations
  #scale_color_manual(values = c('black', 'black', 'black','black','black','black',
   #                            'black','black','black','black','black','steelblue',
    #                           'steelblue','steelblue','steelblue','steelblue','steelblue','steelblue',
     #                          'steelblue','steelblue','steelblue','steelblue','steelblue','steelblue',
      #                         'steelblue','steelblue','steelblue','steelblue','steelblue','black','black','black','black', 'black','black'),
 #                    name = 'Genomic position (3R) 
#Focal allele range: 
#18521012 - 18521208') +
  #geom_count(size = 3) + #Add for points in graph figure if required
  geom_line(linewidth = 1) +
  # Line below can be used to highlight specific datapoint
  #geom_line(data = population %>% filter(Position == '18521012'), color = 'firebrick1') +
  theme_grey() +
  #scale_color_gradient2(low = 'black', mid = 'yellow', high = 'blue') #TODO...gradient code check (relates to line regarding highlighting focal allele)
  labs(title = 'Minor allele frequency within focal region (18520500-18521500) - Cage 1', y = 'Minor allele frequency', x = 'Generation') + # Change Cage number as required
  scale_x_continuous(breaks = c(2,4,8,12,20,28,36,44,56)) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20, face = 'bold'), 
        plot.title = element_text(size = 25, face = 'bold'), 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 20, face = 'bold'))
