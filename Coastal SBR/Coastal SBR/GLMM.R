# Coastal SBR Germination Analysis Workflow

rm(list=ls())
## Load Packages

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(car)
library(emmeans)
library(performance)
library(patchwork)
library(see)


## Import Data

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


## Basic Exploration


summary(dat$GermProportion)

hist(
  dat$GermProportion,
  breaks = 20,
  main = "Germination Proportion",
  xlab = "Proportion Germinated")

mean(dat$GermProportion == 0)


## Germination by Seed Source


ggplot(dat,aes(SeedSource,GermProportion)) +
  geom_boxplot() +
  geom_jitter(
    width = 0.15,
    alpha = 0.4) +
  theme_bw()


## Germination by Soil

ggplot( dat,aes(Soil,GermProportion)) +
  geom_boxplot() +
  geom_jitter(
    width = 0.15,
    alpha = 0.4) +
  theme_bw() +
  coord_flip()


## Germination by Watering


ggplot( dat,aes(Watering, GermProportion)) +
  geom_boxplot() +
  geom_jitter(
    width = 0.15,
    alpha = 0.4) +
  theme_bw()

## Germination by Microtopography


ggplot(dat,aes( Microtopography,GermProportion)) +
  geom_boxplot() +
  geom_jitter(
    width = 0.15,
    alpha = 0.4) +
  theme_bw()

## Continuous Covariates


ggplot(dat,aes( SoilOrganicMatter, GermProportion )) +
  geom_point() +
  geom_smooth(method = "loess" ) +
  theme_bw()


ggplot(dat,aes(SoilBulkDensity,GermProportion)) +
  geom_point() +
  geom_smooth(method = "loess") +
  theme_bw()

ggplot(dat,aes(SoilVWaterContent,GermProportion)) +
  geom_point() +
  geom_smooth( method = "loess") +
  theme_bw()

ggplot(dat,aes(Elevation,GermProportion)) +
  geom_point() +
  geom_smooth( method = "loess") +
  theme_bw()


## Summary Table


dat %>%
  group_by() SeedSource, Soil, Watering, Microtopography) %>%
  summarise(
    n = n(),
    Mean = mean(GermProportion),
    SD = sd(GermProportion),
    Median = median(GermProportion) )


#### Exp 1: Summer Interns ####

exp1<-dat %>% filter(Experiment == "Summer23")

fullmodel1<- glmmTMB(cbind(SeedsGerminated, SeedsSown - SeedsGerminated)~ Phase + Soil , 
                     family = betabinomial(link="logit"), data = exp1)


fullmodelint1<- glmmTMB(cbind(SeedsGerminated, SeedsSown - SeedsGerminated)~ Phase * Soil , 
                        family = betabinomial(link="logit"), data = exp1)


AIC(fullmodel1, fullmodelint1) ### first model without int is lower

anova(fullmodel1, fullmodelint1)

summary(fullmodel1)
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


#### Exp 2: Winter 24 ####

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



#### Exp 3: Winter 26 - Transplant ####

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


#### Exp 4: Winter 26 - Amendment ####

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


#### Graphs ####

library(ggplot2)
library(emmeans)
library(dplyr)

theme_sbr <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.position = "right"
  )

emm_phase1 <- emmeans(fullmodel1, ~ Phase, type = "response")

phase1_df <- as.data.frame(emm_phase1)

gg_exp1 <- ggplot(exp1,
                  aes(x = Phase,
                      y = GermProportion)) +
  geom_jitter(width = 0.15,
              alpha = 0.5,
              size = 2) +
  
  geom_point(
    data = phase1_df,
    aes(x = Phase, y = prob),
    inherit.aes = FALSE,
    size = 4,
    color = "black"
  ) +
  
  geom_errorbar(
    data = phase1_df,
    aes(x = Phase,
        ymin = asymp.LCL,
        ymax = asymp.UCL),
    inherit.aes = FALSE,
    width = 0.15
  ) +
  
  labs(
    x = "Restoration Phase",
    y = "Germination Proportion",
    title = "Experiment 1"
  ) +
  theme_sbr

gg_exp1


emm_soil2 <- emmeans(model2_soil, ~ Soil, type = "response")
soil2_df <- as.data.frame(emm_soil2)

gg_exp2 <- ggplot(exp2,
                  aes(x = Soil,
                      y = GermProportion)) +
  geom_jitter(width = 0.15,
              alpha = 0.5,
              size = 2) +
  
  geom_point(
    data = soil2_df,
    aes(x = Soil,
        y = prob),
    inherit.aes = FALSE,
    size = 4,
    color = "black"
  ) +

  
  labs(
    x = "Soil Type",
    y = "Germination Proportion",
    title = "Experiment 2"
  ) +
  theme_sbr +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))

gg_exp2


emm_source3 <- emmeans(
  fullmodel3,
  ~ SeedSource,
  type = "response"
)

source3_df <- as.data.frame(emm_source3)

gg_exp3 <- ggplot(exp3,
                  aes(x = SeedSource,
                      y = GermProportion)) +
  geom_jitter(width = 0.15,
              alpha = 0.5,
              size = 2) +
  
  geom_point(
    data = source3_df,
    aes(x = SeedSource,
        y = prob),
    inherit.aes = FALSE,
    size = 4,
    color = "black"
  ) +
  
  geom_errorbar(
    data = source3_df,
    aes(x = SeedSource,
        ymin = asymp.LCL,
        ymax = asymp.UCL),
    inherit.aes = FALSE,
    width = 0.15
  ) +
  
  labs(
    x = "Seed Source",
    y = "Germination Proportion",
    title = "Experiment 3"
  ) +
  facet_grid(~Site)+
  theme_sbr

gg_exp3
