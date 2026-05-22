#### Set wd ####

rm(list = ls())
setwd("~/Hester/Coastal SBR/Coastal SBR")


#### load library ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(emmeans)

#### pull in file ####
march<-read.csv('MarchDataSheet.csv')

march <- march %>% 
  mutate(
  SeedSource = factor(SeedSource),
  SiteSown = factor(SiteSown),
  SeedCount = as.numeric(SeedCount))

dataready <- march %>%
  mutate(
    successes = SeedCount,
    failures  = 134 - successes
  ) %>%
  filter(!is.na(successes), SeedCount > 0, successes >= 0, failures >= 0)

#### Early summaries ####


data_summary <- dataready %>%
  mutate(germ = successes) %>%
  group_by(SiteSown, SeedSource) %>%
  summarise(
    n        = n(),
    mean_germ = mean(germ, na.rm = TRUE),
    sd_germ   = sd(germ, na.rm = TRUE),
    se_germ   = sd_germ / sqrt(n)) %>%
  ungroup()

print(data_summary)




data_summary_prop <- dataready %>%
  mutate(germ_prop = successes / 134) %>%
  group_by(SiteSown,SeedSource) %>%
  summarise(
    n        = n(),
    mean_germ = mean(germ_prop, na.rm = TRUE),
    sd_germ   = sd(germ_prop, na.rm = TRUE),
    se_germ   = sd_germ / sqrt(n)
  ) %>%
  ungroup()


Summary(data_summary_prop)


#### make a graph rq ####

ggplot(data_summary_prop,
       aes(x = SeedSource,y = mean_germ, fill = SeedSource)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7,color = "black") +
  geom_errorbar(aes(ymin = mean_germ - se_germ,
                    ymax = mean_germ + se_germ),
                width = 0.2, position = position_dodge(width = 0.8)) +
  facet_wrap(~SiteSown) +
  labs(
    x = "Seed Source",y = "Mean Germination Proportion"  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "none")

