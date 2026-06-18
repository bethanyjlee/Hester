# Coastal SBR Germination Analysis Workflow

## Load Packages

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(car)
library(emmeans)
library(performance)
library(patchwork)
```

## Import Data

```r
dat <- read.csv("CompleteSheet.csv")

dat <- dat %>%
  mutate(
    GermProportion = SeedsGerminated / SeedsSown,
    Soil = factor(Soil),
    Containment = factor(Containment),
    Watering = factor(Watering),
    Microtopography = factor(Microtopography),
    SeedSource = factor(SeedSource),
    Site = factor(Site),
    Experiment = factor(Experiment),
    Phase = factor(Phase)
  )


## Basic Exploration


summary(dat$GermProportion)

hist(
  dat$GermProportion,
  breaks = 20,
  main = "Germination Proportion",
  xlab = "Proportion Germinated"
)

mean(dat$GermProportion == 0)


## Germination by Seed Source


ggplot(
  dat,
  aes(
    SeedSource,
    GermProportion
  )
) +
  geom_boxplot() +
  geom_jitter(
    width = 0.15,
    alpha = 0.4
  ) +
  theme_bw()


## Germination by Soil

ggplot(
  dat,
  aes(
    Soil,
    GermProportion
  )
) +
  geom_boxplot() +
  geom_jitter(
    width = 0.15,
    alpha = 0.4
  ) +
  theme_bw() +
  coord_flip()
```

## Germination by Watering


ggplot(
  dat,
  aes(
    Watering,
    GermProportion
  )
) +
  geom_boxplot() +
  geom_jitter(
    width = 0.15,
    alpha = 0.4
  ) +
  theme_bw()
```

## Germination by Microtopography


ggplot(
  dat,
  aes(
    Microtopography,
    GermProportion
  )
) +
  geom_boxplot() +
  geom_jitter(
    width = 0.15,
    alpha = 0.4
  ) +
  theme_bw()
```

## Continuous Covariates


ggplot(dat,
       aes(
         SoilOrganicMatter,
         GermProportion
       )) +
  geom_point() +
  geom_smooth(
    method = "loess"
  ) +
  theme_bw()


ggplot(dat,
       aes(
         SoilBulkDensity,
         GermProportion
       )) +
  geom_point() +
  geom_smooth(
    method = "loess"
  ) +
  theme_bw()

ggplot(dat,
       aes(
         SoilVWaterContent,
         GermProportion
       )) +
  geom_point() +
  geom_smooth(
    method = "loess"
  ) +
  theme_bw()

ggplot(dat,
       aes(
         Elevation,
         GermProportion
       )) +
  geom_point() +
  geom_smooth(
    method = "loess"
  ) +
  theme_bw()
```

## Summary Table


dat %>%
  group_by(
    SeedSource, Soil, Watering, Microtopography
  ) %>%
  summarise(
    n = n(),
    Mean = mean(GermProportion),
    SD = sd(GermProportion),
    Median = median(GermProportion)
  )
```

## Beta-Binomial Model

```r
m_bb <- glmmTMB(
  cbind(
    SeedsGerminated,
    SeedsSown - SeedsGerminated
  ) ~
    SeedSource +
    Soil +
    Watering +
    Microtopography +
    Containment,
  family = betabinomial(),
  data = dat
)

summary(m_bb)
```

## Type III Tests


car::Anova(
  m_bb,
  type = 3
)
```

## Diagnostics

sim <- simulateResiduals(
  m_bb,
  plot = TRUE
)

testDispersion(sim)

testZeroInflation(sim)

testUniformity(sim)
```

## Estimated Marginal Means

emmeans(m_bb,
  pairwise ~ SeedSource,
  type = "response")

emmeans(m_bb,
  pairwise ~ Soil,
  type = "response")

emmeans(
  m_bb,
  pairwise ~ Watering,
  type = "response"
)
```

## Model Performance

```r
performance::check_model(m_bb)
```

```r
performance::r2(m_bb)
```

