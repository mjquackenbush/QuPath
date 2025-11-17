library(dplyr)
library(ggplot2)
library(stringr)

df <- read_csv(file = "bcellagg_results.csv", locale=locale(encoding="latin1"))
# Filter for only rows with BcellAggregate
df_bcell <- df %>%
  filter(Name == "BcellAggregate")

# Extract week number from Image column (e.g., WK1 -> 1)
df_bcell <- df_bcell %>%
  mutate(
    Week = str_extract(Image, "WK\\d+") %>% str_remove("WK") %>% as.integer()
  )

# Summarize total B cell aggregate area per week
area_summary <- df_bcell %>%
  group_by(Week) %>%
  summarise(TotalArea = sum(`Area µm^2`, na.rm = TRUE))

# Plot
ggplot(area_summary, aes(x = Week, y = TotalArea)) +
  geom_line(group = 1, color = "blue") +
  geom_point(size = 3, color = "darkblue") +
  labs(
    title = "Total Area of B Cell Aggregates Over Time",
    x = "Week",
    y = expression("Total Area ("*µm^2*")")
  ) +
  theme_minimal()
