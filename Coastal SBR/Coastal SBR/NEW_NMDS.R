library(tidyverse)
library(vegan)
library(cluster)
library(dplyr)

rm(list = ls())

setwd("~/GitHub/Hester/Coastal SBR/Coastal SBR")

#### Insert Data ####

dat <- read.csv("CompleteSheet.csv")

#### check variables ####

str(dat)

env <- dat %>%
  select(
    Soil,
    Containment,
    Watering,
    Microtopography,
    SeedSource,
    SoilOrganicMatter,
    SoilBulkDensity,
    SoilVWaterContent,
    Elevation,
    MonthPrecipitation)

dat$GermProportion <- as.numeric(
  as.character(dat$GermProportion))

env <- env %>%
  mutate(
    Soil = factor(Soil),
    Containment = factor(Containment),
    Watering = factor(Watering),
    Microtopography = factor(Microtopography),
    SeedSource = factor(SeedSource))



#### Gowers Conversion ####


gower_dist <- daisy(
  env,
  metric = "gower")



#### starting time ####

nmds <- metaMDS(
  gower_dist,
  k = 2,
  trymax = 500)


cat("Stress =", nmds$stress, "\n")

adonis2(
  gower_dist ~ Soil +
    Watering +
    Microtopography +
    SeedSource +
    Containment,
  data = env
)





scores_df <- as.data.frame(scores(nmds))

plot(
  scores_df$NMDS1,
  scores_df$NMDS2,
  pch = 19,
  col = as.factor(dat$Soil),
  xlab = "NMDS1",
  ylab = "NMDS2"
)

legend(
  "bottomright",
  legend = levels(as.factor(dat$Soil)),
  col = 1:length(levels(as.factor(dat$Soil))),
  pch = 19
)





###########################
# FIT GERMINATION
###########################

fit <- envfit(
  nmds,
  dat["GermProportion"],
  permutations = 999
)

fit

plot(
  fit,
  col = "red"
)

###########################
# FIT ENVIRONMENTAL VARIABLES
###########################

fit_all <- envfit(
  nmds,
  env,
  permutations = 999,
  na.rm = TRUE)

fit_all

plot(nmds, type = "n")
points(scores(nmds), pch = 16)

plot(
  fit_all,
  display = "vectors",
  add = TRUE,
  col = "blue",
  cex = 1,
  arrow.mul = 2
)





#### try again ####
fit_cont <- envfit(
  nmds,
  env %>%
    select(
      Elevation,
      SoilOrganicMatter,
      SoilBulkDensity,
      SoilVWaterContent
    ),
  permutations = 999,
  na.rm = TRUE
)

plot(nmds, type = "n")
points(scores(nmds), pch = 16)

plot(
  fit_cont,
  add = TRUE,
  col = "blue",
  cex = 1.1)
