library(ggplot2)
library(dplyr)

# Folder containing your 9 CSV files
folder <- "/home/baron/Documents/PhD/Data/pop_size_analysis/cage_2/ac_files/final"

# Get all CSV files in that folder
files <- list.files(
  path = folder,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Custom labels for legend in desired order
custom_labels <- c("Gen2", "Gen4", "Gen8", "Gen12", "Gen20",
                   "Gen28", "Gen36", "Gen44", "Gen56")

# Function to read + bin a single file
process_file <- function(f, bin_size = 100) {
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
  geom_line(alpha = 0.8) +
  labs(title = "Binned Line Plot of 9 CSV Files (100 rows per bin)",
       x = "Bin (100 rows each)",
       y = "Mean of Second Column",
       color = "Generation") +  # Updated legend title
  theme_minimal()
