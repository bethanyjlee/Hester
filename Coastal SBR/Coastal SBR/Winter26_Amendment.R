rm(list=ls())

library(dplyr)
library(glmmTMB)
library(car)
library(broom)
library(emmeans)

#### Exp 4: Winter 26 - Amendment ####

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


exp4<-dat %>% filter(Experiment == "Winter26Amendment")

ambient <- subset(exp4, Soil == "Ambient")
phase3 <- subset(exp4, Phase == "Phase3")

model4Phase<-glm(cbind(SeedsGerminated, SeedsSown - SeedsGerminated) ~ Phase,
                 family = quasibinomial, data = ambient)


model4Water<-glm( cbind(SeedsGerminated, SeedsSown - SeedsGerminated) ~ Watering,
                  family = quasibinomial, data = phase3)

Anova(model4Phase, type = "II")
Anova(model4Water, type = "II")

summary(model4Phase)
summary(model4Water)

## Estimated Marginal Means

emmeans(model4Phase,
        pairwise ~ Phase,
        type = "response")

emmeans(model4Water,
        pairwise ~ Watering,
        type = "response")


cont.vars <- exp4 %>%
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
  family = betabinomial(link = "logit"), data = exp4)

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
