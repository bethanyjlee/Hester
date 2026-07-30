# Figure 5A–B: Survival and Growth Barplots with Tukey Letters


rm(list = ls())

setwd("~/Hester/Seed Source/Exp 3 - SoilsGrowth")

library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(survival)
library(survminer)
library(emmeans)
library(multcomp)
library(multcompView)
library(patchwork)

#### READ DATA ####

greenhouse <- read.csv(
  "Elkhorn Greenhouse Data.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE)

names(greenhouse) <- c(
  "TagID", "Soil.Type", "TidalCat", "SiteID",
  "Date1", "Survivors", "Height1",
  "Date2", "Survivors2", "Height2",
  "Date3", "Survivors3", "Height3")

greenhouse$Soil.Type <- recode(
  greenhouse$Soil.Type,
  "Hester" = "Restoration Site")

#### COLORS ####

soil_colors <- c(
  "Restoration Site" = "#999333",
  "Combo" = "#aa4499",
  "Potting Soil" = "#661100")

#### CREATE DAT2 ####

seedsource <- read.csv(
  "FullGermElkhornSeedSource.csv",
  stringsAsFactors = FALSE)

seedsource$TagID <- as.numeric(seedsource$TagID)
greenhouse$TagID <- as.numeric(greenhouse$TagID)

dat2 <- left_join(greenhouse, seedsource, by = "TagID")

#### SURVIVAL SUMMARY ####

surv_plot_data <- dat2 %>%
  mutate(
    Start = as.numeric(Survivors),
    Final = as.numeric(Survivors3),
    Final = pmin(Final, Start),
    SurvivalProp = Final / Start,
    SoilType = factor(
      Soil.Type.x,
      levels = c("Restoration Site", "Combo", "Potting Soil")),
    Tidal = factor(
      TidalCat.x,
      levels = c("Natural", "Diked", "Restored"))) %>%
  filter(
    !is.na(Start),
    !is.na(Final),
    Start > 0)

#### SURVIVAL TUKEY ####

surv_glm <- glm(
  SurvivalProp ~ SoilType * Tidal,
  family = quasibinomial,
  data = surv_plot_data)

emm_surv <- emmeans(surv_glm, ~ SoilType | Tidal)

cld_surv <- cld(emm_surv,Letters = letters,adjust = "tukey")

cld_surv <- as.data.frame(cld_surv)

#### SURVIVAL SUMMARY ####

survival_summary <- surv_plot_data %>%
  group_by(Tidal, SoilType) %>%
  summarise(
    mean_survival = mean(SurvivalProp, na.rm = TRUE),
    se_survival =
      sd(SurvivalProp, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

survival_summary$letters <- cld_surv$.group

#### GROWTH DATA ####

growth_long <- greenhouse %>%
  pivot_longer(
    cols = Date1:Height3,
    names_to = c(".value", "Timepoint"),
    names_pattern = "(Date|Survivors|Height)([0-9]+)") %>%
  mutate( Height = as.numeric(Height),
          SoilType = factor(
            Soil.Type,
            levels = c("Restoration Site", "Combo", "Potting Soil")),
    Tidal = factor(
      TidalCat,
      levels = c("Natural", "Diked", "Restored")))

#### GROWTH MODEL ####

growth_model <- lm(
  Height ~ SoilType * Tidal,
  data = growth_long %>%
    filter(Timepoint == "3"))

emm_growth <- emmeans(growth_model, ~ SoilType | Tidal)

cld_growth <- cld(emm_growth,Letters = letters,adjust = "tukey")

cld_growth <- as.data.frame(cld_growth)

#### GROWTH SUMMARY ####

growth_summary <- growth_long %>%
  filter(Timepoint == "3") %>%
  group_by(Tidal, SoilType) %>%
  summarise(
    mean_growth = mean(Height, na.rm = TRUE),
    se_growth =
      sd(Height, na.rm = TRUE) /
      sqrt(sum(!is.na(Height))),
    .groups = "drop")

growth_summary$letters <- cld_growth$.group

#### FIGURE 5A ####

fig5a <- ggplot(
  survival_summary,
  aes(x = SoilType,y = mean_survival,fill = SoilType)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes
                ( ymin = mean_survival - se_survival,ymax = mean_survival + se_survival),
    width = 0.2,linewidth = 0.8) +
  geom_text(
    aes(
      y = mean_survival + se_survival,
      label = letters
    ),
    vjust = -0.4,
    size = 5
  ) +
  facet_wrap(~ Tidal, ncol = 3,
             labeller = labeller(
               Tidal = c(
                 "Diked" = "Diked Provenances",
                 "Restored" = "Restored Provenances",
                 "Natural" = "Natural Provenances") ) ) + 
  scale_fill_manual(values = soil_colors) +
  ylim(0, 1.05) +
  labs(x = "Soil Treatment",y = "Survival Probability") +
  theme_bw() +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        text = element_text(size = 13))

#### FIGURE 5B ####

fig5b <- ggplot(
  growth_summary,aes
  (x = SoilType,y = mean_growth,fill = SoilType)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes
                (ymin = mean_growth - se_growth, ymax = mean_growth + se_growth),
    width = 0.2,linewidth = 0.8) +
  geom_text(
    aes(
      y = mean_growth + se_growth,
      label = letters),
    vjust = -0.4,
    size = 5) +
  facet_wrap(~ Tidal, ncol = 3,
             labeller = labeller(
               Tidal = c(
                 "Diked" = "Diked Provenances",
                 "Restored" = "Restored Provenances",
                 "Natural" = "Natural Provenances") ) ) + 
  scale_fill_manual(values = soil_colors) +
  labs(x = "Soil Treatment",y = "Final Plant Height (mm)") +
  theme_bw() +
  theme(legend.position = "none",panel.grid.minor = element_blank(),
    text = element_text(size = 13))+
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.10)))

#### combine figure ####

library(cowplot)

final_fig5 <- plot_grid(
  fig5a,
  fig5b,
  labels = c("A", "B"),
  ncol = 1,
  align = "v")

final_fig5
