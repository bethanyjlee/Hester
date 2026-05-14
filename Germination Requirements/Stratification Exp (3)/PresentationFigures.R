# Load the necessary libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(cowplot)  
library(GerminaR)# Load cowplot for easy legend handling

# Step 1: Load the data
new <- read.csv('seedSoaking.csv')

# Step 2: Clean Data - Replace zeros with NA in germination columns
new_clean <- new %>%
  mutate(across(starts_with("D"), ~ ifelse(. == 0, NA, .)))

# Step 3: Recalculate summaries for Potting Soil and Hester Mimic soil
PSgsm_clean <- ger_summary(SeedN = "seeds", evalName = "D", data = new_clean %>% filter(Soil.Type == "Potting Soil"))
HMgsm_clean <- ger_summary(SeedN = "seeds", evalName = "D", data = new_clean %>% filter(Soil.Type == "Hester Mimic"))

# Step 4: Ensure `TimeF` and `SoakF` are factors with correct levels (including Control)
PSgsm_clean$TimeF <- as.factor(PSgsm_clean$Soak.Time)
PSgsm_clean$SoakF <- factor(PSgsm_clean$Salinity.Soak, levels = c("Freshwater", "15ppt", "30ppt", "Control"))

HMgsm_clean$TimeF <- as.factor(HMgsm_clean$Soak.Time)
HMgsm_clean$SoakF <- factor(HMgsm_clean$Salinity.Soak, levels = c("Freshwater", "15ppt", "30ppt", "Control"))

# Step 5: Set the black and white color palette for the plots
bw_palette <- c("black", "grey20", "grey60", "grey80")  # Changed "Control" to grey

# Step 6: Plot 1 - Mean Germination Time for Potting Soil
PSMGT <- ggplot(data = PSgsm_clean, aes(x = TimeF, y = mgt, fill = SoakF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) + 
  scale_fill_manual(values = bw_palette) +
  labs(title = "Days to Germinate in High Nutrient, Potting Soil", 
       x = "Soaking Time", 
       y = "Days to Germinate", 
       fill = "Soaking Salinity") +  # Updated legend title
  scale_y_continuous(limits = c(0, 8)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw() +
  theme(legend.position = "none")  # Remove legend from this plot

# Step 7: Plot 2 - Mean Germination Time for Hester Mimic Soil
HMMGT <- ggplot(data = HMgsm_clean, aes(x = TimeF, y = mgt, fill = SoakF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) + 
  scale_fill_manual(values = bw_palette) +
  labs(title = "Days to Germinate in Nutrient Poor, Hester Soil", 
       x = "Soaking Time", 
       y = "Days to Germinate", 
       fill = "Soaking Salinity") +  # Updated legend title
  scale_y_continuous(limits = c(0, 8)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw() +
  theme(legend.position = "none")  # Remove legend from this plot

# Step 8: Plot 3 - Mean Germinated Seeds for Potting Soil
PSgrc <- ggplot(data = PSgsm_clean, aes(x = TimeF, y = grs, fill = SoakF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) + 
  scale_fill_manual(values = bw_palette) +
  labs(title = "Seeds Germinated in High Nutrient, Potting Soil", 
       x = "Soaking Time", 
       y = "Number of Germinated Seeds", 
       fill = "Soaking Salinity") +  # Updated legend title
  scale_y_continuous(limits = c(0, 35)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw() +
  theme(legend.position = "none")  # Remove legend from this plot

# Step 9: Plot 4 - Mean Germinated Seeds for Hester Mimic Soil
HMgrs <- ggplot(data = HMgsm_clean, aes(x = TimeF, y = grs, fill = SoakF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) + 
  scale_fill_manual(values = bw_palette) +
  labs(title = "Seeds Germinated in Nutrient Poor, Hester Soil", 
       x = "Soaking Time", 
       y = "Number of Germinated Seeds", 
       fill = "Soaking Salinity") +  # Updated legend title
  scale_y_continuous(limits = c(0, 35)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw() +
  theme(legend.position = "none")  # Remove legend from this plot

# Step 10: Create a unified legend
# Create a dummy plot to get the legend
legend <- get_legend(
  ggplot(data = PSgsm_clean, aes(x = TimeF, y = mgt, fill = SoakF)) +
    stat_summary(fun = "mean", geom = "bar", position = position_dodge(width = 0.8), width = 0.7) + 
    scale_fill_manual(values = bw_palette) +
    labs(fill = "Soaking Salinity") +  # Updated legend title
    theme(legend.position = "right")
)

# Step 11: Combine all the plots and place the legend on the right side
combined_plot <- plot_grid(
  PSMGT, HMMGT, PSgrc, HMgrs, 
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
