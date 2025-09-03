# Load libraries
library(ggplot2)
library(readr)
library(dplyr)
library(RColorBrewer)  # for color palettes

# Read in your file
df <- read_table("/home/baron/Documents/PhD/Data/pop_size_analysis/cage_2/AF_calculation/temp/1.txt")  

# Make sure numeric
df$position <- as.numeric(df$position)
df$minor_allele <- as.numeric(df$minor_allele)

# Order chromosomes
chrom_order <- unique(df$chromosome)
df$chromosome <- factor(df$chromosome, levels = chrom_order)

# Compute cumulative positions
df_cum <- df %>%
  group_by(chromosome) %>%
  summarise(chr_len = max(position)) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0)) 

df <- df %>%
  left_join(df_cum, by = "chromosome") %>%
  mutate(pos_cum = position + chr_start)

# Compute chromosome centers for x-axis labels
axis_df <- df %>%
  group_by(chromosome) %>%
  summarise(center = (min(pos_cum) + max(pos_cum)) / 2)

# Generate a distinct color for each chromosome
chrom_colors <- brewer.pal(n = nlevels(df$chromosome), name = "Set1")

# Plot
ggplot(df, aes(x = pos_cum, y = minor_allele, color = chromosome)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_y_continuous(limits = c(0, 0.5)) +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$chromosome
  ) +
  scale_color_manual(values = chrom_colors) +
  labs(
    x = "Chromosome",
    y = "Minor allele frequency"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )