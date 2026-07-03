rm(list=ls())

#### Exp 3: Winter 26 - Transplant ####

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


exp3<-dat %>% filter(Experiment == "Winter26Transplant")

fullmodel3<- glmmTMB(cbind(SeedsGerminated, SeedsSown - SeedsGerminated)~ Site + SeedSource, 
                     family = betabinomial(link="logit"), data = exp3)

fullmodelint3<- glmmTMB(cbind(SeedsGerminated, SeedsSown - SeedsGerminated)~ Site * SeedSource, 
                        family = betabinomial(link="logit"), data = exp3)

AIC(fullmodel3, fullmodelint3) ## better with interaction

anova(fullmodelint3) ## doesnt work so no interactin

Anova(fullmodel3, type = "II")
summary(fullmodel3)

## Diagnostics

sim <- simulateResiduals(
  fullmodel3,
  plot = TRUE)

testDispersion(sim)

testZeroInflation(sim)

testUniformity(sim)


## Estimated Marginal Means

emmeans(fullmodel3,
        pairwise ~ Site,
        type = "response")

emmeans(fullmodel3,
        pairwise ~ SeedSource,
        type = "response")
