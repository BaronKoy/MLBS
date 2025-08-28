# Module load
library(readr)
library(ggplot2)
library(dplyr)

options(scipen = 999) # prevent scientific notation for genomic position on x axis

population <- read_csv('/home/baron/Documents/PhD/Data/plot_data/cage_1/cage_1.csv', # Change as required
                   col_types = cols(Position = col_character(),
                                    Allele_frequency = col_number(),
                                    Generation = col_number()))

# Convert Position to numeric so we can compare ranges
population <- population %>%
  mutate(Position_num = as.numeric(Position),
         Highlight = ifelse(Position_num >= 18521012 & Position_num <= 18521208,
                            "focal", "outside"))

# Plot
ggplot(population, aes(x = Generation, y = Allele_frequency, 
                       group = Position, color = Highlight)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("outside" = "black", "focal" = "steelblue"),
                     name = "Genomic range",
                     labels = c("18521012–18521208 (Suspected BS)", "Positions outside of suspected BS")) +
  labs(title = 'Minor allele frequency within focal region (18520500–18521500) - Cage 1',
       y = 'Minor allele frequency', x = 'Generation') +
  scale_x_continuous(breaks = c(2,4,8,12,20,28,36,44,56)) +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20, face = 'bold'), 
        plot.title = element_text(size = 25, face = 'bold'), 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 20, face = 'bold'))