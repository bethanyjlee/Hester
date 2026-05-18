### Categorical Analysis 


rm(list =ls())
#### set wd

setwd("~/R data/Chp2SeedSource/Exp 2 - Moisture")


### load library ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(emmeans)


###  Read + reshape (cumulative counts) ####
fldata <- read.csv("MoistureStress.csv")

fllong <- fldata %>%
  pivot_longer(cols = starts_with("D"),
               names_to = "Days",
               values_to = "value") %>%
  mutate(
    Days = as.numeric(gsub("D", "", Days)),
    Watering = factor(Watering),
    Seed.Source = factor(Seed.Source),
    Tidal = factor(TidalCat),
    SiteID = factor(SiteID))

dataready <- fllong %>%
  mutate(
    successes = value,
    failures  = Seeds - value
  ) %>%
  filter(!is.na(successes), Seeds > 0, successes >= 0, failures >= 0)

#### Starting from fllong, compute cumulative germination per replicate
fllong_cumgerm <- fllong %>%
  arrange(IDTag, Days) %>%
  group_by(IDTag, SiteID, Watering, Seeds, Tidal) %>%
  mutate(
    alive_count = value,
    cum_germ    = cummax(alive_count),   # running max = total germinated so far
    failures    = Seeds - cum_germ
  ) %>%
  ungroup() %>%
  filter(!is.na(cum_germ), Seeds > 0, cum_germ >= 0, failures >= 0)


#### general summary data ####

data_summary_prop <- dataready %>%
  mutate(germ_prop = successes / Seeds) %>%
  group_by(SiteID, Watering) %>%
  summarise(
    n        = n(),
    mean_germ = mean(germ_prop, na.rm = TRUE),
    sd_germ   = sd(germ_prop, na.rm = TRUE),
    se_germ   = sd_germ / sqrt(n)
  ) %>%
  ungroup()

print(data_summary_prop)


#### select final total germination only for models ####

final_day <- dataready %>%
  group_by(IDTag, SiteID, Watering, TidalCat) %>%
  slice_max(order_by = Days, n = 1, with_ties = FALSE) %>%
  ungroup()

### model assessment and determination ####


mod_final_rand <- glmmTMB(cbind(successes, failures) ~ Watering + TidalCat + (1|SiteID),
                     data = final_day, family = betabinomial())

summary(mod_final_rand)

mod_final_rand_water <- glmmTMB(cbind(successes, failures) ~ Watering + (1|SiteID),
                            data = final_day, family = betabinomial())

summary(mod_final_rand_water)


mod_final_nosouce_int <- glmmTMB(cbind(successes, failures) ~ Watering * Tidal,
                            data = final_day, family = betabinomial())

summary(mod_final_nosouce_int)

mod_final_nosouce_noint <- glmmTMB(cbind(successes, failures) ~ Watering + Tidal,
                                 data = final_day, family = betabinomial())

summary(mod_final_nosouce_noint)

mod_final_sourcet <- glmmTMB(cbind(successes, failures) ~ SiteID,
                              data = final_day, family = betabinomial())

summary(mod_final_sourcet)

####final model assessment - use this one!####

mod_final_int <- glmmTMB(cbind(successes, failures) ~  Watering * TidalCat + (1|SiteID),
                         data = final_day, family = betabinomial())

summary(mod_final_int)


mod_site_water <- glmmTMB(
  cbind(successes, failures) ~ Watering * SiteID,
  data = final_day,
  family = betabinomial())

summary(mod_site_water)


#### Post hoc analyses ####
emm_site_water <- emmeans(
  mod_site_water,
  ~ Watering | SiteID,
  type = "response")

summary(emm_site_water)

pairs(emm_site_water, adjust = "tukey")

emm_sites_by_water <- emmeans(
  mod_site_water,
  ~ SiteID | Watering,
  type = "response")

pairs(emm_sites_by_water, adjust = "tukey")

emm_tidal <- emmeans(
  mod_final_int,
  ~ Watering * TidalCat,
  type = "response")

summary(emm_tidal)

pairs(emm_tidal, adjust = "tukey")

#### tentative plots ####

group_colors<- c(Tidal = "#1E88E5", Muted = '#D81B60', Restored = '#FFC107')
watering_colors <- c (Low = 'orange', High = 'skyblue')

neworder<-c("Tidal","Muted","Restored")
fllong_cumgerm<-arrange(transform(fllong_cumgerm, Tidal = factor(Tidal, levels = neworder)), Tidal)


plot_cum_total_tidal <- fllong_cumgerm %>%
  mutate(prop = cum_germ / Seeds) %>%
  group_by(Watering, SiteID, Days, Tidal) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    se_prop   = sd(prop, na.rm = TRUE) / sqrt(sum(!is.na(prop))),
    .groups   = "drop") %>%
  mutate(lwr = pmax(mean_prop - 1.96 * se_prop, 0),
    upr = pmin(mean_prop + 1.96 * se_prop, 1))




p_spaghetti_flip_total_water <- ggplot() +
 # geom_point(data = plot_cum_total_tidal,
 #            aes(x = Days, y = mean_prop, color = Watering),alpha = 0.7, size = 2)+
  geom_smooth(data = plot_cum_total_tidal, 
              aes(x=Days, y=mean_prop, fill = Watering, color = Watering))+
  facet_wrap(~ Tidal, ncol = 3) +
  scale_y_continuous(limits = c(0, 0.25)) +
  scale_fill_manual(values = watering_colors)+
  scale_color_manual(values = watering_colors) +
  labs(x = "Day",
       y = "Proportion Germinated",
       title = "Cumulative total germination by watering treatment and tidal site history") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

p_spaghetti_flip_total_water

ggsave(filename = "Figure3.tif", dpi = 300, path = "Figures")



