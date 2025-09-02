# Load libraries
library(ggplot2)
library(dplyr)
library(readr)

# Read the file
df <- read_table("/home/baron/Documents/PhD/Data/pop_size_analysis/cage_2/AF_calculation/temp/2.txt")

# Make sure columns are correct
colnames(df) <- c("chromosome", "position", "minor_allele")

# Order chromosomes (important for plotting)
chrom_order <- c("2L", "2R", "3L", "3R", "X")
df$chromosome <- factor(df$chromosome, levels = chrom_order)

# Compute cumulative position across chromosomes
df <- df %>%
  group_by(chromosome) %>%
  mutate(chr_len = max(position)) %>%
  ungroup() %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  group_by(chromosome) %>%
  mutate(pos_cum = position + tot) %>%
  ungroup()

# ---- BINNING (100 rows per bin) ----
df_binned <- df %>%
  group_by(chromosome) %>%
  arrange(pos_cum) %>%
  mutate(bin = (row_number() - 1) %/% 500) %>%   # 100 rows per bin
  group_by(chromosome, bin) %>%
  summarise(
    pos_cum = mean(pos_cum),
    minor_allele = mean(minor_allele),
    .groups = "drop"
  )

# Create axis labels at the midpoint of each chromosome
axis_df <- df_binned %>% 
  group_by(chromosome) %>% 
  summarise(center = (max(pos_cum) + min(pos_cum)) / 2)

# Plot (using the binned data)
ggplot(df_binned, aes(x = pos_cum, y = minor_allele, color = chromosome)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_x_continuous(label = axis_df$chromosome, breaks = axis_df$center) +
  scale_y_continuous(name = "Minor Allele Frequency") +
  labs(x = "Chromosome") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )