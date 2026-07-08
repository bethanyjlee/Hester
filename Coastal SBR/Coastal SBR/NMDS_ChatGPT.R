library(tidyverse)
library(vegan)
library(cluster)

dat <- read.csv("CompleteSheet.csv")

###########################
# CHECK VARIABLES
###########################

str(dat)

###########################
# CLEAN RESPONSE VARIABLE
###########################

dat$GermProportion <- as.numeric(
  as.character(dat$GermProportion)
)

###########################
# SELECT VARIABLES
###########################

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
    MonthPrecipitation
  )

###########################
# CONVERT CATEGORICALS
###########################

env <- env %>%
  mutate(
    Soil = factor(Soil),
    Containment = factor(Containment),
    Watering = factor(Watering),
    Microtopography = factor(Microtopography),
    SeedSource = factor(SeedSource)
  )

###########################
# REMOVE BAD ROWS
###########################

keep <- complete.cases(env) &
  !is.na(dat$GermProportion)

env <- env[keep, ]
dat <- dat[keep, ]

env <- dat %>%
  select(
    Soil,
    Containment,
    Watering,
    Microtopography,
    SeedSource
  )

cat("Rows remaining:", nrow(env), "\n")

###########################
# CHECK FOR PROBLEMS
###########################

sum(is.na(env))
sum(is.na(dat$GermProportion))

env <- env %>%
  mutate(across(where(is.character), as.factor))

###########################
# GOWER DISTANCE
###########################

gower_dist <- daisy(
  env,
  metric = "gower"
)

###########################
# NMDS
###########################

nmds <- metaMDS(
  gower_dist,
  k = 2,
  trymax = 500
)

cat("Stress =", nmds$stress, "\n")

###########################
# PLOT NMDS
###########################

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
  "topright",
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
  permutations = 999
)

fit_all

plot(
  fit_all,
  col = "blue"
)

