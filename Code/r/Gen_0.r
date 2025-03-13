#Modules
library(readr)
library(ggplot2)
library(reshape2)
library(gghighlight)
library(magrittr)
library(dplyr) 

options(scipen = 999) #prevent scientific notation for genomic position on x axis

population <- read_csv('/home/baron/Documents/PhD/ML_balancingselection_data/plot_data/Gen_0.csv',
                   col_types = cols(Position = col_number(), Allele_frequency = col_number(),
                                    DGRP_line = col_character()))

#View file if required
#View(population)

#General plot for all variants within focal region
ggplot(population, aes(x = Position, y = Allele_frequency, color = DGRP_line)) +
  scale_color_manual(values = c('steelblue', 'orange'), 
                     name = 'DGRP line') +
  #geom_count(size = 3) + #Add for points in graph figure if required
  geom_line(linewidth = 1) +
  #Highlight fixed region
  geom_point(data = population %>% filter(Position == '18521012'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521014'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521042'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521044'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521055'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521063'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521069'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521085'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521096'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521103'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521122'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521153'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521162'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521173'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521181'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521182'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521196'), color = 'firebrick1', size = 3) +
  geom_point(data = population %>% filter(Position == '18521208'), color = 'firebrick1', size = 3) +
  #Add annotation to start and end of fixed region
  geom_text(data = population %>% filter(Position == '18521012'), label = '18521012', vjust = -1, color = 'black') +
  geom_text(data = population %>% filter(Position == '18521208'), label = '18521208', vjust = -1, color = 'black') +
  theme_grey() +
  labs(title = 'Alternate allele frequency in original DGRP line samples (14 FB & 14 MB)', y = 'Allele frequency', x = 'Genomic position (Chr3R:18520500-18521500)') +
  #scale_x_discrete(breaks = c(18520500, 18520700, 18520900, 18521100, 18521300, 18521500)) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20, face = 'bold'), 
        plot.title = element_text(size = 25, face = 'bold'), 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 20, face = 'bold'))