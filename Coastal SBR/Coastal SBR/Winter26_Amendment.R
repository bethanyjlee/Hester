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


#### make graph####

library(dplyr)
library(ggplot2)
library(scales)

## Summarize Phase data (Ambient soil only)
phase_plot <- ambient %>%
  group_by(Phase) %>%
  summarise(
    Germination = sum(SeedsGerminated) / sum(SeedsSown),
    .groups = "drop"
  ) %>%
  mutate(Treatment = Phase,
         Variable = "Phase")

## Summarize Watering data (Phase 3 only)
water_plot <- phase3 %>%
  group_by(Watering) %>%
  summarise(
    Germination = sum(SeedsGerminated) / sum(SeedsSown),
    .groups = "drop"
  ) %>%
  mutate(Treatment = Watering,
         Variable = "Watering")

plot_dat <- bind_rows(phase_plot, water_plot)

ggplot(plot_dat,
       aes(x = Treatment, y = Germination)) +
  geom_col(fill = "darkseagreen3",
           color = "black",
           width = 0.7) +
  facet_wrap(~Variable, scales = "free_x") +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Observed Germination (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )


#### Individual Soil Variable Analyses ####

## Soil Bulk Density
model_bulkdensity <- glmmTMB(
  cbind(SeedsGerminated, SeedsSown - SeedsGerminated) ~
    SoilBulkDensity,
  family = betabinomial(link = "logit"),
  data = exp4
)

Anova(model_bulkdensity, type = "II")
summary(model_bulkdensity)

## Estimated marginal trend
emtrends(model_bulkdensity,
         ~1,
         var = "SoilBulkDensity")


## Soil Volumetric Water Content
model_watercontent <- glmmTMB(
  cbind(SeedsGerminated, SeedsSown - SeedsGerminated) ~
    SoilVWaterContent,
  family = betabinomial(link = "logit"),
  data = exp4
)

Anova(model_watercontent, type = "II")
summary(model_watercontent)

emtrends(model_watercontent,
         ~1,
         var = "SoilVWaterContent")


## Soil Organic Matter
model_organicmatter <- glmmTMB(
  cbind(SeedsGerminated, SeedsSown - SeedsGerminated) ~
    SoilOrganicMatter,
  family = betabinomial(link = "logit"),
  data = exp4
)

Anova(model_organicmatter, type = "II")
summary(model_organicmatter)

emtrends(model_organicmatter,
         ~1,
         var = "SoilOrganicMatter")
