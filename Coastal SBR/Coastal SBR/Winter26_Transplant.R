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


#### soil data ####

cont.vars <- exp3 %>%
  dplyr::select(
    SoilOrganicMatter,
    SoilBulkDensity,
    SoilVWaterContent,
    Elevation
  )

cor(cont.vars,
    use = "complete.obs",
    method = "pearson")

library(corrplot)

corrplot(cor(cont.vars,
             use = "complete.obs"),
         method = "number")

summary(aov(SoilOrganicMatter ~ Site, data = exp3))
summary(aov(SoilBulkDensity ~ Site, data = exp3))
summary(aov(SoilVWaterContent ~ Site, data = exp3))
summary(aov(Elevation ~ Site, data = exp3))

summary(aov(SoilOrganicMatter ~ SeedSource, data = exp3))

TukeyHSD(aov(SoilOrganicMatter ~ Site, data = exp3))

TukeyHSD(aov(SoilBulkDensity ~ Site, data = exp3))

TukeyHSD(aov(SoilVWaterContent ~ Site, data = exp3))

TukeyHSD(aov(Elevation ~ Site, data = exp3))

exp3 %>%
  group_by(Site) %>%
  summarise(
    SOM = mean(SoilOrganicMatter),
    SOM_sd = sd(SoilOrganicMatter),
    BulkDensity = mean(SoilBulkDensity),
    BulkDensity_sd = sd(SoilBulkDensity),
    VWC = mean(SoilVWaterContent),
    VWC_sd = sd(SoilVWaterContent),
    Elevation = mean(Elevation),
    Elevation_sd = sd(Elevation)
  )


#### makea da graph ####

library(tidyverse)

soil_plot <- tibble(
  Site = c("Estrada","Hester","Shark Flats"),
  SOM = c(51.2, 2.33, 16.3),
  SOM_sd = c(13.7, 0.402, 2.70),
  BulkDensity = c(0.236, 1.52, 0.807),
  BulkDensity_sd = c(0.192, 0.0474, 0.136),
  VWC = c(0.492, 0.248, 0.403),
  VWC_sd = c(0.0449, 0.0442, 0.101)
)

plot_dat <- bind_rows(
  soil_plot %>%
    transmute(
      Site,
      Variable = "Soil organic matter (%)",
      Mean = SOM,
      SD = SOM_sd
    ),
  
  soil_plot %>%
    transmute(
      Site,
      Variable = "Bulk density (g cm⁻³)",
      Mean = BulkDensity,
      SD = BulkDensity_sd
    ),
  
  soil_plot %>%
    transmute(
      Site,
      Variable = "Volumetric water content",
      Mean = VWC,
      SD = VWC_sd
    )
)

ggplot(plot_dat,
       aes(x = Site,
           y = Mean,
           fill = Site)) +
  geom_col(width = 0.7,
           color = "black") +
  geom_errorbar(aes(ymin = Mean - SD,
                    ymax = Mean + SD),
                width = 0.15) +
  facet_wrap(~Variable,
             scales = "free_y",
             ncol = 2) +
  scale_fill_manual(values = c(
    "Estrada" = "#4daf4a",
    "Hester" = "#e41a1c",
    "Shark Flats" = "#377eb8"
  )) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )


library(dplyr)

germ_plot <- exp3 %>%
  group_by(Site) %>%
  summarise(
    Germination = mean(SeedsGerminated / SeedsSown * 100),
    SD = sd(SeedsGerminated / SeedsSown * 100),
    .groups = "drop"
  )

ggplot(germ_plot,
       aes(x = Site,
           y = Germination,
           fill = Site)) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = Germination - SD,
                    ymax = Germination + SD),
                width = 0.15) +
  scale_fill_manual(values = c(
    "Estrada" = "#4daf4a",
    "Hester" = "#e41a1c",
    "Shark Flats" = "#377eb8"
  )) +
  labs(
    x = NULL,
    y = "Germination (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

germ_source <- exp3 %>%
  group_by(SeedSource) %>%
  summarise(
    Germination = mean(SeedsGerminated / SeedsSown * 100),
    SD = sd(SeedsGerminated / SeedsSown * 100),
    .groups = "drop"
  )

ggplot(germ_source,
       aes(x = SeedSource,
           y = Germination,
           fill = SeedSource)) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = Germination - SD,
                    ymax = Germination + SD),
                width = 0.15) +
  labs(
    x = NULL,
    y = "Germination (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

library(dplyr)
library(ggplot2)

germ_plot <- exp3 %>%
  group_by(Site, SeedSource) %>%
  summarise(
    Mean = mean(SeedsGerminated / SeedsSown * 100),
    SD = sd(SeedsGerminated / SeedsSown * 100),
    n = n(),
    SE = SD/sqrt(n),
    .groups = "drop"
  )

ggplot(germ_plot,
       aes(x = Site,
           y = Mean,
           fill = SeedSource)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7,
           color = "black") +
  geom_errorbar(aes(ymin = Mean - SE,
                    ymax = Mean + SE),
                position = position_dodge(width = 0.8),
                width = 0.2) +
  scale_fill_manual(values = c(
    "Estrada" = "#4daf4a",
    "Hester" = "#e41a1c",
    "SharkFlat" = "#377eb8"
  )) +
  labs(
    x = "Transplant Site",
    y = "Germination (%)",
    fill = "Seed Source"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )
