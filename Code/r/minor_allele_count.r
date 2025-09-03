library(ggplot2)
library(dplyr)

# Folder containing your 9 CSV files
folder <- "/home/baron/Documents/PhD/Data/pop_size_analysis/cage_2/ac_files/final"

# Get all CSV files in that folder
files <- list.files(
  path = folder,
  pattern = "\\.txt$", # Change for csv=empirical, txt=simulated - in /home/baron/Documents/PhD/Data/pop_size_analysis/cage_2/ac_files/final
  full.names = TRUE
)

# Custom labels for legend in desired order
custom_labels <- c("Gen2", "Gen56")

# Function to read + bin a single file
process_file <- function(f, bin_size = 500) {
  df <- read.table(f, header = FALSE, sep = "", stringsAsFactors = FALSE)
  colnames(df) <- c("Col1", "Col2")
  
  df <- df %>%
    mutate(Row = seq_along(Col2),
           Value = Col2,
           Bin = Row %/% bin_size) %>%
    group_by(Bin) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(File = basename(f))
  
  return(df)
}

# Process all files and combine
all_data <- bind_rows(lapply(files, process_file))

# Map filenames to custom labels
label_map <- setNames(custom_labels, basename(files))
all_data$FileLabel <- label_map[all_data$File]

# Turn FileLabel into a factor with levels in the order you want
all_data$FileLabel <- factor(all_data$FileLabel, levels = custom_labels)

# Plot with custom legend in numerical order
ggplot(all_data, aes(x = Bin, y = Value, color = FileLabel)) +
  geom_line(alpha = 0.8, linewidth = 2) +
  labs(title = "Mean minor allele count (500 sites per bin) - simulated(generation 2 vs 56)", x = "Bin count (500 sites per bin)",
       y = expression("Mean minor allele count (500 sites)" ~ "=" ~ f[500]),
       color = "Generation") +
  scale_color_manual(
    values = c("Gen2" = "steelblue", "Gen56" = "orange")   # pick your colors here
  ) +
  scale_y_continuous(limits = c(15, 35)) +
  theme_grey() +
  theme(
    plot.title = element_text(size = 30, face = "bold"),       # Title size
    axis.title = element_text(size = 20),                      # Axis labels size
    axis.text = element_text(size = 15),                       # Axis tick labels size
    legend.title = element_text(size = 15),                    # Legend title size
    legend.text = element_text(size = 15)                      # Legend item text size
  )