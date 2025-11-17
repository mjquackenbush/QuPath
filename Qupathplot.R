#QuPath quantification plotting of measurement tables
library(dplyr)
library(ggpubr)
df <- read_tsv(file = "measurement_CD11b+Ly6G+PDL1+.tsv")
# Assuming your dataframe is called df
# Get all column names starting with "Num"
num_cols <- grep("^Num", names(df), value = TRUE)

# Loop through each "Num" column and create new columns
for (col in num_cols) {
  density_col <- paste0("Density_", col)
  proportion_col <- paste0("Proportion_", col)
  
  df[[density_col]] <- df[[col]] / df$Area
  df[[proportion_col]] <- df[[col]] / df$Detections
}


df <- df %>%
  mutate(Condition = case_when(
    grepl("WT", Image) ~ "WT",
    grepl("CD8 ?KO", Image) ~ "CD8 KO",
    grepl("CD8 ST", Image) ~ "CD8 ST",
    grepl("CD8 LT", Image) ~ "CD8 LT",
    TRUE ~ "Other"
  ))

ggplot(df, aes(x = Condition, y = `Proportion_Num CD11b`, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") +
  stat_compare_means(method = "kruskal.test", label.y = max(df$`Proportion_Num CD11b`) * 1.1) +  # global test
  stat_compare_means(
    method = "wilcox.test", 
    comparisons = combn(unique(df$Condition), 2, simplify = FALSE),
    label = "p.signif"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Proportion of CD11b+ Cells by Condition",
    x = "Condition",
    y = "Proportion CD11b"
  ) +
  scale_fill_brewer(palette = "Set2")

#for-loop for all cell combinations in df
# Identify all columns that start with "Proportion_"
proportion_cols <- grep("^Proportion_", names(df), value = TRUE)

# Loop through each proportion column and save plot
for (col in proportion_cols) {
  # Generate the plot title (e.g., "Proportion CD11b")
  # Extract marker name (e.g., from "Proportion_Num CD11b" → "CD11b")
  marker <- gsub("^Proportion_Num\\s*", "", col)
  marker <- gsub("_", " ", marker)  # clean any remaining underscores
  
  # Create the formatted title
  plot_title <- paste("Proportion of", paste0(marker, "+ Cells by Condition"))
  
  
  # Create the plot
  p <- ggplot(df, aes(x = Condition, y = .data[[col]], fill = Condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black") +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") +
    stat_compare_means(method = "kruskal.test", label.y = max(df[[col]], na.rm = TRUE) * 1.1) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = combn(unique(df$Condition), 2, simplify = FALSE),
      label = "p.signif"
    ) +
    theme_minimal(base_size = 14) +
    labs(
      title = plot_title,
      x = "Condition",
      y = plot_title
    ) +
    scale_fill_brewer(palette = "Set2")
  
  # Save the plot as a TIFF
  ggsave(
    filename = paste0(col, "_by_condition.tiff"),
    plot = p,
    device = "tiff",
    width = 7,
    height = 5,
    dpi = 300
  )
}

#Density
library(ggplot2)
library(ggpubr)

# Identify all columns that start with "Density_"
density_cols <- grep("^Density_", names(df), value = TRUE)

# Loop through each density column and save plot
for (col in density_cols) {
  # Extract marker name (e.g., from "Density_Num CD11b" → "CD11b")
  marker <- gsub("^Density_Num\\s*", "", col)
  marker <- gsub("_", " ", marker)  # replace underscores with spaces
  
  # Create the formatted title
  plot_title <- paste("Density of", paste0(marker, "+ Cells by Condition"))
  
  # Create the plot
  p <- ggplot(df, aes(x = Condition, y = .data[[col]], fill = Condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black") +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") +
    stat_compare_means(method = "kruskal.test", label.y = max(df[[col]], na.rm = TRUE) * 1.1) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = combn(unique(df$Condition), 2, simplify = FALSE),
      label = "p.signif"
    ) +
    theme_minimal(base_size = 14) +
    labs(
      title = plot_title,
      x = "Condition",
      y = plot_title
    ) +
    scale_fill_brewer(palette = "Set2")
  
  # Save the plot as a TIFF
  ggsave(
    filename = paste0(col, "_by_condition.tiff"),
    plot = p,
    device = "tiff",
    width = 7,
    height = 5,
    dpi = 300
  )
}


df <- df %>%
  mutate(Region = case_when(
    grepl("TUMOR", Classification, ignore.case = TRUE) ~ "Tumor",
    grepl("Stroma", Classification, ignore.case = TRUE) ~ "Stroma",
    TRUE ~ "Other"
  ))


df$Group <- interaction(df$Condition, df$Region, sep = "_")

ggplot(df, aes(x = Group, y = `Proportion_Num CD11b`, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") +
  stat_compare_means(
    method = "kruskal.test", 
    label.y = max(df$`Proportion_Num CD11b`) * 1.1
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Proportion of CD11b+ Cells by Condition and Region",
    x = "Condition + Region",
    y = "Proportion CD11b"
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = function(x) gsub("_", "\n", x))  # optional: splits label into 2 lines


ggplot(df, aes(x = Region, y = `Proportion_Num CD11b`, fill = Region)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") +
  stat_compare_means(
    aes(group = Region), 
    method = "wilcox.test", 
    label = "p.signif"
  ) +
  facet_wrap(~ Condition, nrow = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Proportion of CD11b+ Cells by Region Within Each Condition",
    x = "Region",
    y = "Proportion CD11b"
  ) +
  scale_fill_brewer(palette = "Set2")


