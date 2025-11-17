#classify cells as near, intermediate, or far based on their distance to B cells
# Load libraries
library(readr)
library(dplyr)
library(stringr)
# Read the exported QuPath file (adjust path/sep if needed)
data <- read_delim("measurementss.txt", delim = "\t")

# Extract Condition and Tumor replicate info
data <- data %>%
  mutate(
    Condition = sub("_.*", "", Image),
    Tumor = str_extract(Image, "T\\d+")
  )

# Filter to CD8a cells with valid distance to B cells
cd8_cells <- data %>%
  filter(Classification == "CD8" & !is.na(`Distance to B cells (µm)`))

# Filter to CD8+ cells with valid distance to B cells
cd8_cells <- data %>%
  filter(Classification == "CD8", !is.na(`Distance to B cells (µm)`)) %>%
  mutate(
    Distance_Class = case_when(
      `Distance to B cells (µm)` <= 10 ~ "Adjacent",
      `Distance to B cells (µm)` <= 50 ~ "Intermediate",
      TRUE ~ "Far"
    ),
    Distance_Class = factor(Distance_Class, levels = c("Adjacent", "Intermediate", "Far"))
  )

# Calculate per-tumor mean distance
tumor_means <- cd8_cells %>%
  group_by(Condition, Tumor) %>%
  summarise(Mean_Distance = mean(`Distance to B cells (µm)`), .groups = "drop")

#plot class distribution
# Plot: boxplot by condition + tumor replicate means
ggplot(cd8_cells, aes(x = Condition, y = `Distance to B cells (µm)`, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_point(data = tumor_means,
             aes(x = Condition, y = Mean_Distance),
             color = "black", size = 3, shape = 21, fill = "white", stroke = 1) +
  labs(title = "CD8+ T Cell Distance to B Cells by Condition",
       x = "Condition",
       y = "Distance to B cells (µm)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")
library(ggpubr)
ggplot(cd8_cells, aes(x = Condition, y = `Distance to B cells (µm)`, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_point(data = tumor_means,
             aes(x = Condition, y = Mean_Distance),
             color = "black", size = 3, shape = 21, fill = "white", stroke = 1) +
  labs(title = "CD8+ T Cell Distance to B Cells by Condition",
       x = "Condition",
       y = "Distance to B cells (µm)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  stat_compare_means(method = "wilcox.test", label.y = max(cd8_cells$`Distance to B cells (µm)`) * 1.05)

library(ggplot2)
# Plot distribution of distance classes by condition
ggplot(cd8_cells, aes(x = Distance_Class, fill = Condition)) +
  geom_bar(position = "dodge") +
  labs(title = "CD8 T cells by distance to nearest B cell",
       x = "Distance Class",
       y = "Cell Count") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")

# Bar plot of CD8+ T cell distance categories
ggplot(cd8_cells, aes(x = Distance_Class)) +
  geom_bar(fill = "#4682B4") +
  labs(title = "CD8 T cells by distance to nearest B cell", x = "Distance Class", y = "Cell Count") +
  theme_minimal()

library(spatstat.geom)
library(spatstat.explore)

# Create spatial point pattern (ppp) for CD8+ T cells
ppp_cd8 <- ppp(
  x = cd8_cells$`Centroid X µm`,
  y = cd8_cells$`Centroid Y µm`,
  window = owin(range(cd8_cells$`Centroid X µm`), range(cd8_cells$`Centroid Y µm`))
)

# Ripley's K function
K <- Kest(ppp_cd8)
plot(K)


# Clustering vs. Randomness
#If your empirical K(r) (any of the colored lines) is above the blue dashed Poisson line, it indicates spatial clustering — more neighbors than expected under randomness.
#If below, it would suggest spatial inhibition — points repel each other.


#Plateauing
#If K(r) flattens out (plateaus), it may indicate that most clustering occurs within a certain distance range.
#Your lines don’t fully plateau, so clustering persists out to the max distance (80 µm), but may be strongest around 10–40 µm.


#Ripleys CrossK
# Load packages
library(spatstat)
library(spatstat.geom)
library(dplyr)

# Read your QuPath export (assumed cleaned)
df <- read.delim("measurementss.txt")

# Filter only Cd8 and Cd20 cells (adjust class column name if needed)
cells <- df %>%
  filter(Class %in% c("Cd8a", "Cd20")) %>%
  mutate(mark = Class)

# Rename coordinate columns
names(cells)[names(cells) == "Centroid.X.µm"] <- "X"
names(cells)[names(cells) == "Centroid.Y.µm"] <- "Y"

# Check counts
table(cells$mark)
# Create point pattern (replace W with correct image size in microns)
# Construct the observation window from min/max centroid positions
library(spatstat)

# Define spatial window in microns
W <- owin(xrange = range(cells$X, na.rm = TRUE),
          yrange = range(cells$Y, na.rm = TRUE))

# Create point pattern object with cell type as marks
pp <- ppp(x = cells$X,
          y = cells$Y,
          window = W,
          marks = as.factor(cells$mark))
# Run Ripley's Cross-K to examine spatial interaction between Cd8a and Cd20
Kcross_result <- Kcross(pp, i = "Cd8a", j = "Cd20", correction = "Ripley")
Lcross_result <- Lcross(pp, i = "Cd8a", j = "Cd20", correction = "Ripley")
# Plot K
plot(Kcross_result, main = "Ripley's Cross-K: CD8 vs B Cells")

# Plot L (interpreted as: L(r) > r means clustering)
plot(Lcross_result, main = "Ripley's L-function: CD8 vs B Cells")
abline(0, 1, col = "red", lty = 2)  # Line of CSR

#L(r) > r: Attraction (more likely to find CD8s near B cells at that radius).
#L(r) < r: Repulsion (CD8s avoid B cells).
#L(r) ≈ r: No interaction (random distribution).


#ask reverse, do B cells cluster around Cd8s?
Lcross_result_BtoT <- Lcross(pp, i = "Cd20", j = "Cd8a", correction = "Ripley")
# Plot L-function for reverse direction
plot(Lcross_result_BtoT, main = "Ripley's L-function: B Cells vs CD8 T Cells")
abline(0, 1, col = "red", lty = 2)  # CSR reference line



#data distribution
ggplot(cd8_cells, aes(x = `Distance to B cells (µm)`, fill = Condition)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 40) +
  facet_wrap(~Condition, scales = "free_y") +
  labs(title = "Histogram of CD8+ T Cell Distances",
       x = "Distance to B cells (µm)",
       y = "Count") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")

#Count CD8+ cells per tumor, per distance class
cd8_counts <- cd8_cells %>%
  group_by(Condition, Tumor, Distance_Class) %>%
  summarise(CD8_Count = n(), .groups = "drop")
#Get total cell counts per tumor (from your full dataset before filtering)
total_cells <- data %>%
  group_by(Condition, Tumor) %>%
  summarise(Total_Cells = n(), .groups = "drop")
#Normalize CD8 counts by total tumor cells
cd8_normalized <- cd8_counts %>%
  left_join(total_cells, by = c("Condition", "Tumor")) %>%
  mutate(CD8_Proportion = CD8_Count / Total_Cells)
#Compare across conditions (Adj, Intermediate)
cd8_subset <- cd8_normalized %>%
  filter(Distance_Class %in% c("Adjacent", "Intermediate","Far"))
#plot
ggplot(cd8_subset, aes(x = Condition, y = CD8_Proportion, fill = Condition)) +
  geom_boxplot(alpha = 0.6) +
  facet_wrap(~Distance_Class, scales = "free_y") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Proportion of CD8+ Cells Adjacent, Intermediate, or Far to B Cells",
       y = "Proportion of CD8+ Cells",
       x = "Condition") +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = list(c("CD40", "NT"))
  )

