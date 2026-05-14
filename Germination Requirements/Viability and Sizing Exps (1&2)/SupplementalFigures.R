library(dplyr)
library(ggplot2)
library(readr)

## -------------------------
## Read your data
## -------------------------
dat <- read_csv("GermSized.csv")

## -------------------------
## Make sure columns are named correctly
## (EDIT if needed)
## -------------------------
# Expected columns:
# SiteID = site identifier (x-axis)
# SeedSize = numeric values (y-axis)

dat <- dat %>%
  rename(
    SiteID = SiteID,      # change if needed
    SeedSize = Size   # change if needed
  )

## -------------------------
## Summary stats for labels
## -------------------------
sum_dat <- dat %>%
  group_by(SiteID) %>%
  summarise(
    mean_size = mean(SeedSize, na.rm = TRUE),
    .groups = "drop"
  )

## OPTIONAL: If you have Tukey letters, merge them here
# letters_df <- read_csv("letters.csv")
# sum_dat <- left_join(sum_dat, letters_df, by = "SiteID")

## -------------------------
## Plot
## -------------------------
p <- ggplot(dat, aes(x = SiteID, y = SeedSize)) +
  geom_boxplot(
    fill = "lightblue",
    color = "black",
    outlier.shape = NA   # removes outlier points
  ) +
  
  # mean labels (add letters if you have them)
  geom_text(
    data = sum_dat,
    aes(
      x = SiteID,
      y = max(dat$SeedSize, na.rm = TRUE) + 0.1,
      label = round(mean_size, 2)
      # if using letters: paste0(round(mean_size, 2), " ", letters)
    ),
    size = 4
  ) +
  
  labs(
    title = "Seed Size Distribution",
    x = "Site ID",
    y = "Seed Size (mm)"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

## -------------------------
## Save figure
## -------------------------
ggsave(
  filename = "LeeSupplementalFigure4.tif",
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)
