rm(list=ls())

library(tidyr)
library(dplyr)
library(glmmTMB)
library(car)
library(broom)
library(emmeans)

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
