# Plot Ne from analysis of single Ne per cage population

# ---- LIBRARY PREPARATION ----
library(tidyverse)
library(ggplot2)

# ---- DATASET CREATION ----
df <- tribble(
  ~Cage_population, ~Ne, ~bootstrap,
  "Cage_1", 253, "248.33 – 253.73 – 258.94",
  "Cage_2", 268, "260.92 – 265.85 – 272.34",
  "Cage_3", 241, "235.71 – 240.95 – 245.58",
  "Cage_4", 219, "217.33 – 222.43 – 227.56",
  "Cage_5", 207, "201.28 – 206.53 – 210.94",
  "Cage_6", 279, "269.06 – 274.22 – 278.52",
  "Cage_7", 297, "288.84 – 295.42 – 300.95",
  "Cage_8", 306, "291.2 – 297.7 – 304.56",
  "Cage_9", 274, "267.5 – 273.55 – 279.97",
  "Cage_10", 304, "292.36 – 298.35 – 303.89"
)

# ---- CLEAN AND PREPARE DATA ----
df_clean <- df %>%
  separate(bootstrap, into = c("lower", "median", "upper"), sep = " – ", convert = TRUE) %>%
  mutate(Cage_population = factor(Cage_population, levels = paste0("Cage_", 1:10)))  # 1 at top

# ---- PLOT ----
ggplot(df_clean, aes(y = Cage_population)) +
  # CI bars
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = "95% Bootstrap CI"), 
                 height = 0.3, linewidth = 1) +
  # Main estimate
  geom_point(aes(x = Ne, color = "Ne estimate"), size = 3.5) +
  # Bootstrap median
  geom_point(aes(x = median, color = "Bootstrap median"), 
             shape = 21, fill = "white", size = 3.5, stroke = 1) +
  
  # Color scheme
  scale_color_manual(
    name = NULL,
    values = c(
      "95% Bootstrap CI" = "gray50",
      "Ne estimate" = "#1F78B4",      # steel blue
      "Bootstrap median" = "#E69F00"  # orange, color-blind safe
    )
  ) +
  
  # Labels
  labs(
    x = expression("Effective Population Size (N"["e"]*")"),
    y = "Cage Population",
    title = "Effective Population Size with Bootstrap Intervals"
  ) +
  
  # Clean theme
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major.x = element_line(color = "gray90", linetype = "dotted"),
    panel.grid.minor.x = element_blank()
  ) +
  coord_flip()