# Load libraries
library(ggplot2)
library(dplyr)
library(readr)
library(purrr)

# ---- SETTINGS ----
# Folder with your .txt files
data_dir <- "/home/baron/Documents/PhD/Data/pop_size_analysis/cage_2/AF_calculation/temp"

# Get all .txt files
files <- list.files(path = data_dir, pattern = "\\.txt$", full.names = TRUE)

# ---- READ + COMBINE ----
all_data <- map_dfr(files, function(f) {
  df <- read_table(f, col_names = TRUE)
  colnames(df) <- c("chromosome", "position", "minor_allele")
  df$file <- basename(f)  # track which file
  df
})

# ---- CHROMOSOME ORDER ----
chrom_order <- c("2L", "2R", "3L", "3R", "X") # adjust if you have more
all_data$chromosome <- factor(all_data$chromosome, levels = chrom_order)

# ---- CUMULATIVE POSITION ----
all_data <- all_data %>%
  group_by(chromosome) %>%
  mutate(chr_len = max(position)) %>%
  ungroup() %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  group_by(chromosome, file) %>%
  mutate(pos_cum = position + tot) %>%
  ungroup()

# ---- BINNING (100 rows per bin) ----
all_binned <- all_data %>%
  group_by(file, chromosome) %>%
  arrange(pos_cum) %>%
  mutate(bin = (row_number() - 1) %/% 500) %>%   # 100 SNPs per bin
  group_by(file, chromosome, bin) %>%
  summarise(
    pos_cum = mean(pos_cum),
    minor_allele = mean(minor_allele),
    .groups = "drop"
  )

# ---- X-AXIS LABELS ----
axis_df <- all_binned %>%
  group_by(chromosome) %>%
  summarise(center = (max(pos_cum) + min(pos_cum)) / 2)

# ---- PLOT WITH CORRECT CHROMOSOME AXIS ----
ggplot(all_binned, aes(x = pos_cum, y = minor_allele, color = file)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(se = FALSE, method = "loess", span = 0.2) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = axis_df$center,              # midpoint per chromosome
    labels = axis_df$chromosome           # show chromosome names
  ) +
  scale_y_continuous(name = "Minor Allele Frequency") +
  labs(color = "File") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )