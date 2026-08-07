rm(list=ls())

library(tidyr)
library(dplyr)
library(glmmTMB)
library(car)
library(broom)
library(emmeans)
library(ggplot2)

#### Exp 1: Summer Interns ####

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
    Phase = factor(Phase))


exp1<-dat %>% filter(Experiment == "Summer23")

fullmodel1<- glmmTMB(cbind(SeedsGerminated, SeedsSown - SeedsGerminated)~ Phase + Soil , 
                     family = betabinomial(link="logit"), data = exp1)


fullmodelint1<- glmmTMB(cbind(SeedsGerminated, SeedsSown - SeedsGerminated)~ Phase * Soil , 
                        family = betabinomial(link="logit"), data = exp1)


AIC(fullmodel1, fullmodelint1) ### first model without int is lower

anova(fullmodel1, fullmodelint1)

summary(fullmodel1)

#### include soil data ####





## Type III Tests

car::Anova(fullmodel1, type = 3)

car::Anova(fullmodel1, type = "II")

## Diagnostics

sim <- simulateResiduals(
  fullmodel1,
  plot = TRUE)

testDispersion(sim)

testZeroInflation(sim)

testUniformity(sim)


## Estimated Marginal Means

emmeans(fullmodel1,
        pairwise ~ Soil,
        type = "response")

emmeans(fullmodel1,
        pairwise ~ Phase,
        type = "response")


## Model Performance


performance::check_model(fullmodel1)

performance::r2(fullmodel1)


#### make graph ####

scale_fill_manual(values = c(
  "Ambient" = "#56B4E9",      # Sky blue
  "Compost" = "#009E73",      # Bluish green
  "PeatPot" = "#E69F00",      # Orange
  "PottingSoil" = "#CC79A7"   # Reddish purple
))

exp1 <- dat %>%
  filter(Experiment == "Summer23")

exp1_summary <- exp1 %>%
  group_by(Soil, Phase) %>%
  summarise(
    mean_germ = mean(GermProportion),
    sd_germ = sd(GermProportion),
    .groups = "drop"
  )

ggplot(exp1_summary, aes(x = Soil, y = mean_germ, fill = Soil)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_germ - sd_germ,
        ymax = mean_germ + sd_germ),
    width = 0.2) +
  geom_jitter(
    data = exp1,
    aes(x = Soil, y = GermProportion),
    inherit.aes = FALSE,
    width = 0.1,
    size = 2) +
  facet_wrap(~Phase)+
  labs(y = "Mean Germination") +
  scale_fill_manual(values = c(
    "Ambient"      = "#56B4E9",
    "Coir"         = "#0072B2",
    "Compost"      = "#009E73",
    "PeatPot"      = "#E69F00",
    "PottingSoil"  = "#CC79A7",
    "YampahMud"    = "#D55E00"))
