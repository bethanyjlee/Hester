rm(list=ls())

library(tidyverse)
library(glmmTMB)
library(car)
library(emmeans)
library(DHARMa)
library(performance)
library(ggplot2)


garden <- read.csv("MoniqueGarden.csv")

# Keep only needed variables
garden <- garden %>%
  select(
    Area,
    Microtopography,
    PercentCover,
    Seedlings
  )

# Clean Seedlings column
garden$Seedlings <- as.character(garden$Seedlings)

garden$Seedlings[garden$Seedlings == "1 large"] <- "1"
garden$Seedlings[garden$Seedlings == "3 large"] <- "3"

garden$Seedlings <- as.numeric(garden$Seedlings)

garden$Area <- factor(garden$Area)
garden$Microtopography <- factor(garden$Microtopography)


# model discovery
m1 <- glmmTMB(
  Seedlings ~ Area * Microtopography,
  family = poisson,
  data = garden
)

simulationOutput <- simulateResiduals(m1)
plot(simulationOutput)

check_overdispersion(m1)

m2 <- glmmTMB(
  Seedlings ~ Area * Microtopography,
  family = nbinom2(),
  data = garden
)

Anova(m2, type = "III")
summary(m2)

simulationOutput <- simulateResiduals(m2)
plot(simulationOutput)

check_overdispersion(m2)

m3 <- update(m2, . ~ Area + Microtopography)

Anova(m3, type = "II")
summary(m3)

AIC(m2, m3)

# emmeans

emm_area <- emmeans(m3, ~ Area)

emm_area

pairs(emm_area)

emm_micro <- emmeans(m3, ~ Microtopography)

emm_micro

pairs(emm_micro)

# plot raw data
emm_plot <- as.data.frame(emmeans(m3, ~ Area, type = "response"))

garden <- garden %>%
  filter(Microtopography != "") %>%
  droplevels()

#plot 1
ggplot(garden, aes(x = Area, y = Seedlings, color = Area)) +
  geom_jitter(
    width = 0.15,
    height = 0,
    size = 3,
    alpha = 0.7) +
  geom_point(
    data = emm_plot,
    aes(x = Area, y = response),
    inherit.aes = FALSE,
    size = 5,
    color = "black") +
  geom_errorbar(
    data = emm_plot,
    aes(x = Area,y = response,
      ymin = asymp.LCL,
      ymax = asymp.UCL ),
    inherit.aes = FALSE,
    width = 0.15,
    linewidth = 1,
    color = "black") +
  scale_color_manual(values = c(
    "Control" = "#0072B2",
    "SeedAddition" = "#D55E00" )) +
  labs( x = "",
    y = "Number of Seedlings") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")

# plot 2
ggplot(garden, aes(Area, Seedlings, fill = Area)) +
  geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.7) +
  scale_fill_manual(values = c(
    "Control" = "#0072B2",
    "SeedAddition" = "#D55E00"
  )) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none") +
  labs(x = "",y = "Seedling Count")

#plot 3
ggplot(garden, aes(Area, Seedlings, fill = Area)) +
  geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.7) +
  scale_fill_manual(values = c(
    "Control" = "#0072B2",
    "SeedAddition" = "#D55E00"
  )) +
  facet_wrap(~Microtopography) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none") +
  labs(x = "",y = "Seedling Count")
