# Load libraries
library(ggplot2)
library(dplyr)
library(readr)
library(purrr)

# ---- SETTINGS ----
data_dir <- "/home/baron/Documents/PhD/Data/pop_size_analysis/cage_2/AF_calculation/temp"

# Get all .txt files
files <- list.files(path = data_dir, pattern = "\\.txt$", full.names = TRUE)

# ---- Pick one random SNP position from the first file ----
first_df <- read_table(files[1], col_names = TRUE)
colnames(first_df) <- c("chromosome", "position", "minor_allele")

random_snp <- first_df %>%
  slice_sample(n = 1) %>%
  select(chromosome, position)

print(random_snp)  # to see which SNP was chosen

# ---- Extract this SNP across all files ----
all_snps <- map_dfr(files, function(f) {
  df <- read_table(f, col_names = TRUE)
  colnames(df) <- c("chromosome", "position", "minor_allele")
  df$file <- basename(f)
  df %>% semi_join(random_snp, by = c("chromosome", "position"))
})

# ---- RENAME + ORDER FILES FOR X-AXIS ----
# Create a lookup table with descriptive names
file_labels <- c(
  "1.txt" = "Gen 2",
  "2.txt" = "Gen 4",
  "3.txt" = "Gen 8",
  "4.txt" = "Gen 12",
  "5.txt" = "Gen 20",
  "6.txt" = "Gen 28",
  "7.txt" = "Gen 36",
  "8.txt" = "Gen 44",
  "9.txt" = "Gen 56"
)

# Apply labels
all_snps$file <- file_labels[all_snps$file]

# Ensure correct order
all_snps$file <- factor(all_snps$file, levels = file_labels)

# ---- PLOT ----
ggplot(all_snps, aes(x = file, y = minor_allele, group = 1, color = file)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1, color = "black", alpha = 0.6) +
  labs(
    x = "Generation",
    y = "Minor Allele Frequency",
    title = paste0("Minor Allele Frequency for SNP ",
                   random_snp$chromosome, ":", random_snp$position)
  ) +
  scale_y_continuous(limits = c(0, 0.50)) +   # <--- fixed scale
  theme_gray() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )