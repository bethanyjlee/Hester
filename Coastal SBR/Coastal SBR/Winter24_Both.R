rm(list=ls())

#### Exp 2: Winter 24 ####

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


exp2<-dat %>% filter(Experiment == "Winter24")

fullmodel2<- glmmTMB(cbind(SeedsGerminated, SeedsSown - SeedsGerminated)~ Soil + Microtopography, 
                     family = betabinomial(link="logit"), data = exp2)

fullmodelint2<- glmmTMB(cbind(SeedsGerminated, SeedsSown - SeedsGerminated)~ Soil * Microtopography, 
                        family = betabinomial(link="logit"), data = exp2)

m_soil2 <- glmmTMB( cbind(SeedsGerminated, SeedsSown - SeedsGerminated) ~ Soil,
                    family = betabinomial(), data = exp2)

m_cont2 <- glmmTMB( cbind(SeedsGerminated, SeedsSown - SeedsGerminated) ~ Containment,
                    family = betabinomial(), data = exp2)

AIC(fullmodel2, m_cont2, m_soil2) ### AIC is better for all but model design was poor so we must analyze separately

anova(fullmodel2, fullmodelint2)

summary(m_soil2)
summary(m_cont2)


### not gonna use betabinomial for this experiment

kruskal.test(GermProportion ~ Soil, data = exp2)

kruskal.test(GermProportion ~ Containment, data = exp2)

model2_soil<-glm(cbind(SeedsGerminated, SeedsSown - SeedsGerminated) ~ Soil,
                 family = quasibinomial, data = exp2)

model2_cont<-glm(cbind(SeedsGerminated, SeedsSown - SeedsGerminated) ~ Containment,
                 family = quasibinomial, data = exp2)

## Type III Tests


car::Anova(model2_cont, type = 3)
car::Anova(model2_soil, type = 3)


## Estimated Marginal Means

emmeans(model2_cont,
        pairwise ~ Containment,
        type = "response")

emmeans(model2_soil,
        pairwise ~ Soil,
        type = "response")



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

