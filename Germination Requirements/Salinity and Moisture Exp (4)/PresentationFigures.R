# Load the necessary libraries
library(ggplot2)
library(GerminaR)
library(dplyr)
library(gridExtra)
library(cowplot)  # for handling the shared legend

# Step 1: Load the data
germsum <- read.csv('Mega2New.csv')

# Step 2: Clean Data - Replace zeros with NA in germination columns
germsum_clean <- germsum %>%
  mutate(across(starts_with("D"), ~ ifelse(. == 0, NA, .)))

# Step 3: Recalculate summaries for Potting Soil and Hester Mimic soil
PSgsm_clean <- ger_summary(SeedN = "Seeds", evalName = "D", data = germsum_clean %>% filter(Soil == "PottingSoil"))
HMgsm_clean <- ger_summary(SeedN = "Seeds", evalName = "D", data = germsum_clean %>% filter(Soil == "HesterMimic"))

# Step 4: Ensure `TimeF` and `SoakF` are factors with correct levels
PSgsm_clean$MoistureF <- factor(PSgsm_clean$Moisture, levels = c("Low", "Medium", "High"))
HMgsm_clean$MoistureF <- factor(HMgsm_clean$Moisture, levels = c("Low", "Medium", "High"))

# Set color palette
bw_palette <- c("black", "grey20", "grey60", "grey80")  # Ensuring "Control" is grey

# Step 5: Plot 1 - Germinated Seeds in Potting Soil (PSGRS)
PSGRS <- ggplot(data = PSgsm_clean, aes(x = Salinity, y = grs, fill = MoistureF)) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = bw_palette) +
  labs(title = "Seeds Germinated in High Nutrient, Potting Soil", x = "Moisture Addition to Soils", 
       y = "Number of Germinated Seeds", fill = "Moisture Levels") +
  scale_y_continuous(limits = c(0, 35)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw() +
  theme(legend.position = "none")  # Remove legend from this plot

# Step 6: Plot 2 - Mean Germination Time in Potting Soil (PSMGT)
PSMGT <- ggplot(data = PSgsm_clean, aes(x = Salinity, y = mgt, fill = MoistureF)) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = bw_palette) +
  labs(title = "Days to Germinate in High Nutrient, Potting Soil", x = "Moisture Addition to Soils", 
       y = "Days to Germinate", fill = "Moisture Levels") +
  scale_y_continuous(limits = c(0, 16)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw() +
  theme(legend.position = "none")  # Remove legend from this plot

# Step 7: Plot 3 - Germinated Seeds in Hester Mimic Soil (HMGRS)
HMGRS <- ggplot(data = HMgsm_clean, aes(x = Salinity, y = grs, fill = MoistureF)) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = bw_palette) +
  labs(title = "Seeds Germinated in Nutrient Poor, Hester Soil", x = "Moisture Addition to Soils", 
       y = "Number of Germinated Seeds", fill = "Moisture Levels") +
  scale_y_continuous(limits = c(0, 35)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw() +
  theme(legend.position = "none")  # Remove legend from this plot

# Step 8: Plot 4 - Mean Germination Time in Hester Mimic Soil (HMMGT)
HMMGT <- ggplot(data = HMgsm_clean, aes(x = Salinity, y = mgt, fill = MoistureF)) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = bw_palette) +
  labs(title = "Days to Germinate in Nutrient Poor, Hester Soil", x = "Moisture Addition to Soils", 
       y = "Days to Germinate", fill = "Moisture Levels") +
  scale_y_continuous(limits = c(0, 16)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw() +
  theme(legend.position = "none")  # Remove legend from this plot

# Step 9: Create a unified legend
# Create a dummy plot to get the legend
legend <- get_legend(
  ggplot(data = PSgsm_clean, aes(x = Salinity, y = mgt, fill = MoistureF)) +
    stat_summary(fun = "mean", geom = "bar", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = bw_palette) +
    labs(fill = "Moisture Levels") +  # Unified legend title
    theme(legend.position = "right")
)

# Step 10: Combine all the plots and place the legend on the right side
combined_plot <- plot_grid(
  PSMGT, HMMGT,   # First row: Germination Time (Potting Soil & Hester Soil)
  PSGRS, HMGRS,   # Second row: Germinated Seeds (Potting Soil & Hester Soil)
  ncol = 2, 
  labels = c("A", "B", "C", "D"), 
  label_size = 12,
  rel_widths = c(1, 1), 
  rel_heights = c(1, 1)
)

# Combine the plots with the legend
final_plot <- plot_grid(
  combined_plot, legend, 
  ncol = 2, 
  rel_widths = c(1, 0.2)  # Adjust the width to ensure the legend fits on the right
)

# Display the final plot
print(final_plot)
