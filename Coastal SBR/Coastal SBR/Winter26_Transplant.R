rm(list=ls())

library(dplyr)
library(glmmTMB)
library(car)
library(broom)
library(emmeans)

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

anova(fullmodelint3) ## doesnt work so no interaction

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

#### Model Discovery with numeric variables ####

## Continuous variables
cont.vars <- exp3 %>%
  dplyr::select(SoilOrganicMatter,
                SoilBulkDensity,
                SoilVWaterContent,
                Elevation,
                MonthPrecipitation)

## Correlation matrix
cor(cont.vars,
    use = "complete.obs",
    method = "pearson")

## Optional visualization
library(corrplot)

corrplot(cor(cont.vars,
             use = "complete.obs"),
         method = "number")


envmodel <- glmmTMB(
  cbind(SeedsGerminated, SeedsSown - SeedsGerminated) ~
    Site + SeedSource + SoilOrganicMatter +
    SoilBulkDensity + SoilVWaterContent + Elevation +MonthPrecipitation,
  family = betabinomial(link = "logit"), data = exp3)

Anova(envmodel, type = "II")

summary(envmodel)

AIC(fullmodel3, envmodel)

### estimated marginal means

emmeans(envmodel,
        pairwise ~ Site,
        type = "response")

emmeans(envmodel,
        pairwise ~ SeedSource,
        type = "response")
 
#### model selection ####
m0 <- glmmTMB(
  cbind(SeedsGerminated, SeedsSown-SeedsGerminated) ~
    Site + SeedSource,
  family=betabinomial(link="logit"),
  data=exp3)

m1 <- update(m0, . ~ . + SoilOrganicMatter)
m2 <- update(m0, . ~ . + SoilBulkDensity)
m3 <- update(m0, . ~ . + SoilVWaterContent)
m4 <- update(m0, . ~ . + Elevation)

AIC(m0, m1, m2, m3, m4)

summary(m4)
