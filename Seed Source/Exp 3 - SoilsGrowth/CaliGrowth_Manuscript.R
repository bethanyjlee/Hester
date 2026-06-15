rm(list = ls())

setwd("~/Hester/Seed Source/Exp 3 - SoilsGrowth")

library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)


#### Read greenhouse data ####
greenhouse <- read.csv(
  "Elkhorn Greenhouse Data.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE)

#### Reshape to LONG format ####

# Step 1: Fix duplicate column names
names(greenhouse) <- c(
  "TagID", "Soil.Type", "TidalCat", "SiteID",
  "Date1", "Survivors", "Height1",
  "Date2", "Survivors2", "Height2",
  "Date3", "Survivors3", "Height3")

# Step 2: Reshape into long format
growth_long <- greenhouse %>%
  pivot_longer(cols = Date1:Height3,
    names_to = c(".value", "Timepoint"), 
    names_pattern = "(Date|Survivors|Height)([0-9]+)") %>%
  mutate(TagID = as.numeric(TagID),
    Height = as.numeric(Height),
    Survivors = as.numeric(Survivors)) %>%
  mutate(SoilType = factor(Soil.Type, levels = c("Hester", "Combo", "Potting Soil")),
    SiteID = factor(SiteID),
    Tidal = factor(TidalCat),
    Timepoint = factor(Timepoint,
                       levels = c("1","2","3"),
                       labels = c("6/12", "7/22", "8/21")),
    Date = factor(Date))


#### Quick Plot: growth by siteID #### 
ggplot(growth_long, aes(x = SiteID, y = Height, fill = SoilType)) +
  stat_summary(fun = mean, geom = "bar",
    position = position_dodge(width = 0.8), width = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.8),width = 0.2) +
  labs(x = "Seed Source",y = "Mean Tallest Plant Height (mm)",
       title = "Plant Growth by Tidal Category") +
  theme_bw()




#### Final Height Summary Table by Seed Source ####

final_growth_summary <- growth_long %>%
  filter(Timepoint == "8/21") %>%
  group_by(SiteID) %>%
  summarise(
    Total_Plants = sum(!is.na(Height)),
    Mean_Final_Height = mean(Height, na.rm = TRUE),
    SD_Final_Height = sd(Height, na.rm = TRUE),
    SE_Final_Height = SD_Final_Height / sqrt(Total_Plants)
  ) %>%
  arrange(SiteID)

print(final_growth_summary)

#### Model Discovery --- SiteID ####

dataready <- growth_long %>%
  mutate(
    Soils = factor(SoilType, levels = c("Hester", "Combo", "Potting Soil")),
    SiteId = factor(SiteID),
    growth = factor(Height),
    Tidal = factor(TidalCat, levels = c("Natural", "Diked", "Restored")))

# Model with interactions
m_full <- glm(
  growth ~ Soils * SiteID,
  data = dataready,
  family = binomial)
summary(m_full)

# Model without interactions
m_no_interaction <- glm(
  growth ~ Soils + SiteID,
  data = dataready,
  family = binomial)


anova(m_full, m_no_interaction, test = "Chisq") # AIC is lower for with interactions
# use m_full but lets check family use


#quasi 
m_quasi <- glm(
  growth ~ Soils * SiteID,
  data = dataready,
  family = quasibinomial)
summary(m_quasi)
car::Anova(m_quasi, type = "II")


#betabinomial
m_bb <- glmmTMB(
  growth ~ Soils * SiteID,
  data = dataready,
  family = betabinomial())
summary(m_bb)

#### Model Discovery --- tidal categories ####

# Model with interactions
m_full_tid <- glm(
  growth ~ Soils * Tidal,
  data = dataready,
  family = binomial)
summary(m_full_tid)

# Model without interactions
m_no_interaction_tid <- glm(
  growth ~ Soils + Tidal,
  data = dataready,
  family = binomial)
summary(m_no_interaction_tid)

anova(m_full_tid, m_no_interaction_tid, test = "Chisq") # AIC is lower for with interactions
# use m_full but lets check family use


#quasi 
m_quasi_tid <- glm(
  growth ~ Soils * Tidal,
  data = dataready,
  family = quasibinomial)

summary(m_quasi_tid)
car::Anova(m_quasi_tid, type = "II")


#betabinomial
m_bb_tid <- glmmTMB(
  growth ~ Soils * Tidal,
  data = dataready,
  family = betabinomial())

summary(m_bb_tid)
car::Anova(m_bb_tid, type = "II")



library(emmeans)

emmeans(m_bb_tid, pairwise ~ Soils, type = "response")

emmeans(m_bb_tid, pairwise ~ Soils | Tidal, type = "response")


#### Post hoc analyses ####


#### OVERALL MODEL TESTS — TIDAL ####

library(car)

# Type III Wald tests
car::Anova(m_bb_tid, type = "III")

# pseudo-R2
library(performance)
r2(m_bb_tid)

# confidence intervals
confint(m_bb_tid)

# likelihood ratio test for interaction
m_bb_tid_no_int <- glmmTMB(
  growth ~ Soils + Tidal,
  data = dataready,
  family = betabinomial()
)

anova(m_bb_tid_no_int, m_bb_tid)

#### EMMEANS — SOIL EFFECTS ####

library(emmeans)

emm_growth_tid <- emmeans(
  m_bb_tid,
  ~ Soils | Tidal,
  type = "response"
)

summary(emm_growth_tid)

# pairwise soil comparisons within tidal category
growth_pairs <- pairs(
  emm_growth_tid,
  by = "Tidal",
  adjust = "tukey"
)

summary(growth_pairs, infer = TRUE)

#### SITE EFFECTS ####

car::Anova(m_bb, type = "III")

# emmeans by site
emm_growth_site <- emmeans(
  m_bb,
  ~ SiteID,
  type = "response"
)

pairs(emm_growth_site, adjust = "tukey")

#### ESTIMATED MEANS ####

summary(emm_growth_tid) %>%
  as.data.frame()




#### full plots ###
group_colors<- c(Natural = "#1E88E5", Diked = '#D81B60', Restored = '#FFC107')
soil_colors<- c(Hester = "#999333", Combo = '#aa4499', "Potting Soil" = '#661100')

neworder<-c("Natural","Diked","Restored")
growth_long<-arrange(transform(growth_long, Tidal = factor(Tidal, levels = neworder)), Tidal)



## Plot: growth by TidalCat ##
ggplot(growth_long, aes(x = Tidal, y = Height, fill = SoilType)) +
  stat_summary(fun = mean, geom = "bar",
               position = position_dodge(width = 0.8), width = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.8),width = 0.2) +
  labs(x = "Seed Source",y = "Mean Tallest Plant Height (mm)",
       title = "Plant Growth by Tidal Category") +
  scale_fill_manual(values = soil_colors) +
  theme_bw()

## Plot: growth by Site Faceted by Soil ##
ggplot(growth_long, aes(x = Tidal, y = Height, fill = Tidal)) +
  stat_summary(fun = mean, geom = "bar",
               position = position_dodge(width = 0.8), width = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.8), width = 0.2) +
  scale_fill_manual(values = group_colors) +   
  labs(x = "Seed Source",
       y = "Mean Tallest Plant Height (mm)",
       title = "Plant Growth by Tidal Category") +
  facet_wrap(~SoilType) +
  theme_bw()

## Plot: growth by Indiviudal Site 
ggplot(growth_long, aes(x = SiteID, y = Height, fill = SoilType)) +
  stat_summary(fun = mean, geom = "bar",
               position = position_dodge(width = 0.8), width = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.8),width = 0.2) +
  labs(x = "Seed Source",y = "Mean Tallest Plant Height (mm)",
       title = "Plant Growth by Tidal Category") +
  scale_fill_manual(values = soil_colors) +
  theme_bw()
