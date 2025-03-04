# Modules
library(readr)
library(ggplot2)
library(reshape2)
library(gghighlight)
options(scipen = 999) # prevent scientific notation for genomic position on x axis

population <- read_csv('/home/baron/Documents/PhD/ML_balancingselection_data/plot_data/region_data_1.csv', #Change file as required
                   col_types = cols(Position = col_character(), Allele_frequency = col_number(),
                                    Generation = col_number()))
#focal <- population$Position['18521001'] #TODO...Code for highlighting FOCAL
#View file if required
#View(allele)

# General plot for all variants within focal region
ggplot(population, aes(x = Generation, y = Allele_frequency, color = Position)) +
  #geom_count(size = 3) + #Add for points in graph figure if required
  geom_line(linewidth = 0.5) +
  #gghighlight(max(Allele_frequency) > 0.3) + #Highlights for threshold frequencies, change value as required
  theme_grey() +
  #scale_color_gradient2(low = 'black', mid = 'yellow', high = 'blue') #TODO...gradient code check (relates to line regarding highlighting focal allele)
  labs(title = 'Alternate allele frequency within focal region (18520500-18521500) - Cage 1', y = 'Alternate allele frequency', x = 'Generation') + # Change Cage number as required
  scale_x_continuous(breaks = c(2,4,8,12,20,28,36,44,56)) +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))
