rm(list=ls())


library(dplyr)
library(ggplot2)
library(tidyverse)

dat<-read.csv("WalledGardens.csv")

long_dat <- dat %>%
  pivot_longer(
    cols = c(SeedlingsQ1, AdultsQ1,
             SeedlingsQ2.5, AdultsQ2.5,
             SeedlingsQ4, AdultsQ4),
    names_to = c("LifeStage", "Quadrat"),
    names_pattern = "(Seedlings|Adults)Q(.+)",
    values_to = "Count"
  ) %>%
  mutate(
    Quadrat = factor(Quadrat,
                     levels = c("1", "2.5", "4")),
    Walled = factor(Walled,
                    levels = c("No", "Yes"))
  )

head(long_dat)

glm(Count ~ Walled, family = poisson(), 
    data = filter(long_dat, LifeStage == "Seedlings"))

library(glmmTMB)

mod<-glmmTMB(
  Count ~ Walled + (1 | `Number..E.to.W.`),
  family = poisson(),
  data = filter(long_dat, LifeStage == "Seedlings")
)

car::Anova(mod)

library(ggplot2)

seedlings <- long_dat %>%
  filter(LifeStage == "Seedlings")

ggplot(seedlings,
       aes(x = Walled, y = Count, fill = Walled)) +
  stat_summary(fun = mean,
               geom = "bar",
               width = 0.6,
               color = "black") +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.15,
               linewidth = 0.7) +
  geom_jitter(width = 0.08,
              size = 2,
              alpha = 0.7) +
  labs(x = "Walled or no?", y = "Mean seedling abundance" ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13))

ggplot(seedlings,
       aes(Walled, Count, color = Walled)) +
  geom_jitter(width = 0.12,
              size = 2.5,
              alpha = 0.8) +
  stat_summary(fun = mean,
               geom = "point",
               size = 5,
               color = "black") +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.12,
               color = "black") +
  labs(x = NULL,y = "Seedling abundance" ) +
  theme_classic() +
  theme(legend.position = "none")

library(dplyr)

seedlings <- long_dat %>%
  filter(LifeStage == "Seedlings")

seedlings %>%
  group_by(Walled) %>%
  summarise(
    mean = mean(Count),
    sd = sd(Count),
    se = sd(Count)/sqrt(n()),
    median = median(Count),
    min = min(Count),
    max = max(Count),
    n = n()
  )


mod.nb <- glmmTMB(
  Count ~ Walled + (1 | Number..E.to.W.),
  family = nbinom2(),
  data = filter(long_dat, LifeStage == "Seedlings")
)

summary(mod.nb)
AIC(mod, mod.nb)






##### Survival of species #####
library(tidyverse)

surv <- dat %>%
  select(
    Garden = `Number..E.to.W.`,
    Walled,
    Distichlis,
    Frankenia,
    Jaumea
  ) %>%
  pivot_longer(
    cols = c(Distichlis, Frankenia, Jaumea),
    names_to = "Species",
    values_to = "Alive"
  ) %>%
  mutate(
    Dead = 5 - Alive,
    Survival = Alive / 5,
    Garden = factor(Garden),
    Walled = factor(Walled, levels = c("No", "Yes")),
    Species = factor(Species)
  )

head(surv)

library(glmmTMB)

mod_surv <- glmmTMB(
  cbind(Alive, Dead) ~ Walled * Species +
    (1 | Garden),
  family = binomial,
  data = surv
)

summary(mod_surv)
car::Anova(mod_surv)

library(ggplot2)

ggplot(surv,
       aes(x = Walled, y = Survival, fill = Walled)) +
  stat_summary(fun = mean,
               geom = "bar",
               width = 0.7,
               color = "black") +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.2) +
  geom_jitter(width = 0.08,
              size = 2,
              alpha = 0.7) +
  facet_wrap(~Species) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Walled or No?",
    y = "Survival proportion"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

mod_surv2 <- glmmTMB(
  cbind(Alive, Dead) ~ Walled + Species +
    (1 | Garden),
  family = binomial,
  data = surv
)

summary(mod_surv2)
car::Anova(mod_surv2)

library(emmeans)

emmeans(mod_surv2, pairwise ~ Species, type = "response")

emmeans(mod_surv2, pairwise ~ Walled, type = "response")

## analyze each species seprately
mod.frank <- glmmTMB(
  cbind(Alive, Dead) ~ Walled + (1|Garden),
  family = binomial,
  data = filter(surv, Species == "Frankenia")
)

car::Anova(mod.frank)

mod.dist <- glmmTMB(
  cbind(Alive, Dead) ~ Walled + (1|Garden),
  family = binomial,
  data = filter(surv, Species == "Distichlis")
)

car::Anova(mod.dist)

mod.jaum <- glmmTMB(
  cbind(Alive, Dead) ~ Walled + (1|Garden),
  family = binomial,
  data = filter(surv, Species == "Jaumea")
)

car::Anova(mod.jaum)



library(tidyr)

paired_seedlings <- seedlings %>%
  select(`Number..E.to.W.`, Quadrat, Walled, Count) %>%
  pivot_wider(
    names_from = Walled,
    values_from = Count
  ) %>%
  drop_na()

# If approximately normal
t.test(
  paired_seedlings$Yes,
  paired_seedlings$No,
  paired = TRUE
)

# If not normal (recommended for count data)
wilcox.test(
  paired_seedlings$Yes,
  paired_seedlings$No,
  paired = TRUE,
  exact = FALSE
)



##### Effect of groundwater influence #####

seedlings <- seedlings %>%
  mutate(
    GroundwaterInfluenced = factor(Groundwater.influenced))

mod.gw <- glmmTMB(
  Count ~ Walled + GroundwaterInfluenced +
    (1 | Number..E.to.W.),
  family = nbinom2(),
  data = seedlings
)

summary(mod.gw)
car::Anova(mod.gw)

mod.gw.int <- glmmTMB(
  Count ~ Walled * GroundwaterInfluenced +
    (1 | Number..E.to.W.),
  family = nbinom2(),
  data = seedlings
)

summary(mod.gw.int)
car::Anova(mod.gw.int)

library(emmeans)
emmeans(mod.gw.int, pairwise ~ Walled | GroundwaterInfluenced,
        type = "response")

with(seedlings, table(GroundwaterInfluenced, Walled))

aggregate(
  Count ~ GroundwaterInfluenced + Walled,
  data = seedlings,
  summary
)

#### okay since we have one corner of our 2x2 treatments that is 0 lets just subset

gw <- subset(seedlings, GroundwaterInfluenced == "Yes")

glmmTMB(
  Count ~ Walled + (1|Number..E.to.W.),
  family = nbinom2(),
  data = gw
)

##### Vegetation (natural seedlings only) by garden

garden_seedlings <- seedlings %>%
  group_by(`Number..E.to.W.`, Walled) %>%
  summarise(
    TotalSeedlings = sum(Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Garden = factor(`Number..E.to.W.`,
                    levels = sort(unique(`Number..E.to.W.`)))
  )

ggplot(garden_seedlings,
       aes(x = Garden,
           y = TotalSeedlings,
           fill = Walled)) +
  geom_col(position = position_dodge(width = 0.8),
           color = "black",
           width = 0.7) +
  labs(
    x = "Garden",
    y = "Naturally Recruited Seedlings",
    fill = "Walled"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title = element_text(size = 13)
  )

##### Vegetation (natural seedlings only) by groundwater-influenced gardens

garden_seedlings <- seedlings %>%
  filter(`Number..E.to.W.` %in% 6:10) %>%
  group_by(`Number..E.to.W.`, Walled) %>%
  summarise(
    TotalSeedlings = sum(Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Garden = factor(`Number..E.to.W.`, levels = 6:10)
  )

ggplot(garden_seedlings,
       aes(x = Garden,
           y = TotalSeedlings,
           fill = Walled)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7,
           color = "black") +
  labs(
    x = "Groundwater-influenced garden",
    y = "Total naturally recruiting seedlings",
    fill = "Walled"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 13),
    legend.position = "top"
  )
