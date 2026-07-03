rm(list=ls())

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

