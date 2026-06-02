# Attempting a PCA for all SBR efforts
# started May 22nd, 2026

library(cluster)
library(tidyverse)

rm(list = ls())

setwd("~/Hester/Coastal SBR/Coastal SBR")

####pull in data ####


dat <- read.csv("IncompleteSheet.csv")

pca_data <- dat %>%
  select(
    SoilOrganicMatter,
    SoilBulkDensity,
    SoilVWaterContent,
    Elevation,
    MonthPrecipitation,
    AnnualPrecipitation,
  )

# Remove rows with missing values
pca_data <- na.omit(pca_data)

# Run PCA
pca <- prcomp(
  pca_data,
  center = TRUE,
  scale. = TRUE
)

summary(pca)

# Scree plot
plot(pca)

# Biplot
biplot(pca)