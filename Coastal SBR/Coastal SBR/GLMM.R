library(tidyverse)
library(lme4)
library(car)
library(DHARMa)
library(performance)
library(MuMIn)
library(emmeans)

rm(list = ls())

setwd("~/GitHub/Hester/Coastal SBR/Coastal SBR")

dat <- read.csv("CompleteSheet.csv")

str(dat)

dat <- dat %>%
  mutate(
    Soil = factor(Soil),
    Containment = factor(Containment),
    Watering = factor(Watering),
    Microtopography = factor(Microtopography),
    SeedSource = factor(SeedSource),
    Site = factor(Site),
    Experiment = factor(Experiment),
    Phase = factor(Phase)
  )

hist(dat$GermProportion)

summary(dat$GermProportion)


m1 <- glm(
  cbind(
    SeedsGerminated,
    SeedsSown - SeedsGerminated
  ) ~
    SeedSource +
    Soil +
    Watering +
    Microtopography +
    Containment,
  family = binomial,
  data = dat
)

summary(m1)


Anova(
  m1,
  type = 3
)

m_glm <- glm(
  cbind(
    SeedsGerminated,
    SeedsSown - SeedsGerminated
  ) ~
    SeedSource +
    Soil +
    Watering +
    Microtopography +
    Containment,
  family = binomial,
  data = dat
)

dat %>%
  group_by(Soil) %>%
  summarise(
    n = n(),
    Germinated = sum(SeedsGerminated),
    Seeds = sum(SeedsSown),
    Prop = Germinated / Seeds
  )


dat %>%
  group_by(Containment) %>%
  summarise(
    Germinated = sum(SeedsGerminated),
    Seeds = sum(SeedsSown),
    Prop = Germinated / sum(SeedsSown)
  )
