

### set library ###
library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(emmeans)

### set wd ###
setwd("~/R data/Chp2SeedSource/Exp 2 - Moisture")



#### Read + reshape (cumulative germ)####
fldata <- read.csv("MoistureStress.csv")

fllong <- fldata %>%
  pivot_longer(cols = starts_with("D"),
               names_to = "Days",
               values_to = "value") %>%
  mutate(Days = as.numeric(gsub("D", "", Days)),
    Watering = factor(Watering),
    Seed.Source = factor(SiteID),
    TidalCat = factor(TidalCat))

dataready <- fllong %>%
  mutate(successes = value,failures  = Seeds - value) %>%
  filter(!is.na(successes), Seeds > 0, successes >= 0, failures >= 0)


data_summary <- dataready %>%
  mutate(germ_prop = successes / Seeds) %>%
 group_by(Seed.Source, TidalCat) %>%
 # group_by(Watering) %>%
  summarise(
    n        = n(),
    mean_germ = mean(germ_prop, na.rm = TRUE),
    sd_germ   = sd(germ_prop, na.rm = TRUE),
    se_germ   = sd_germ / sqrt(n)) %>%
  ungroup()

print(data_summary)

####plotting ####
### prep ###

group_colors<- c(Tidal = "#1E88E5", Muted = '#D81B60', Restored = '#FFC107')

plot_cum <- dataready %>%
  mutate(prop = successes / Seeds) %>%
  group_by(Watering, Seed.Source, Days, TidalCat) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    se_prop   = sd(prop, na.rm = TRUE) / sqrt(sum(!is.na(prop))),
    .groups = "drop") %>%
  mutate(lwr = pmax(mean_prop - 1.96 * se_prop, 0),
    upr = pmin(mean_prop + 1.96 * se_prop, 1))

final_day <- dataready %>%
  group_by(IDTag, Seed.Source, Watering, TidalCat) %>%
  slice_max(order_by = Days, n = 1, with_ties = FALSE) %>%
  ungroup()


### model selection/summary ###
mod_final_int <- glmmTMB(cbind(successes, failures) ~  Watering * TidalCat + (1|Seed.Source),
                         data = final_day, family = betabinomial())

summary(mod_final_int)


### emmeans info and CI ###
emm_final <- as.data.frame(emmeans(mod_final_int, ~ Watering * TidalCat, type = "response"))

## Robust CI + estimate column handling (fixes your lower.CL error)
emm_final <- emm_final %>%
  mutate(
    prob = dplyr::coalesce(
      if ("prob" %in% names(.)) prob else NA_real_,
      if ("response" %in% names(.)) response else NA_real_,
      if ("rate" %in% names(.)) rate else NA_real_),
    lwr = dplyr::coalesce(
      if ("lower.CL" %in% names(.)) lower.CL else NA_real_,
      if ("asymp.LCL" %in% names(.)) asymp.LCL else NA_real_,
      if ("lower.HPD" %in% names(.)) lower.HPD else NA_real_),
    upr = dplyr::coalesce(
      if ("upper.CL" %in% names(.)) upper.CL else NA_real_,
      if ("asymp.UCL" %in% names(.)) asymp.UCL else NA_real_,
      if ("upper.HPD" %in% names(.)) upper.HPD else NA_real_))


#### OPTION 1: Cumulative curves (mean ± 95% CI) by day #### 

p_cum <- ggplot(plot_cum, aes(x = Days, y = mean_prop, color = Watering, group = Watering)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Watering), alpha = 0.18, color = NA) +
  facet_wrap(~ Seed.Source) +
  scale_y_continuous(limits = c(0, 0.3)) +
  scale_fill_manual(values = c("High" = "blue","Low"  = "red"))+
  scale_color_manual(values = c("High"= "blue", "Low"='red'))+
  labs(x = "Day", y = "Cumulative germination (proportion)",
       title = "Cumulative germination curves") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

p_cum


## -------------------------
## OPTION 2: Final day only (recommended for inference)
##   Uses the last available day per IDTag
## -------------------------


## EMMEANS interaction plot (points + CI)
p_final_emm <- ggplot(emm_final, aes(x = TidalCat, y = prob, color = Watering, group = Watering)) +
  geom_point(position = position_dodge(width = 0.35), size = 2) +
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                position = position_dodge(width = 0.35), width = 0.15) +
  geom_line(position = position_dodge(width = 0.35)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  labs(x = "Seed source", y = "Predicted final germination (proportion)",
       title = "Final-day model EMMEANS (beta-binomial)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
p_final_emm

## Bar version
p_final_bar <- ggplot(emm_final, aes(x = TidalCat, y = prob, fill = Watering)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                position = position_dodge(width = 0.75), width = 0.2) +
  scale_y_continuous(limits = c(0, 0.4)) +
  scale_fill_manual(values = c("High" = "blue","Low"  = "red")) +
  labs(x = "Seed source", y = "Predicted final germination (proportion)",
       title = "Final-day model predictions by treatment") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
p_final_bar

## -------------------------
## OPTION 3 (Flipped):
## Two panels = Watering
## Colored by Seed.Source
## -------------------------

rep_curves <- dataready %>%
  mutate(prop = successes / Seeds)

mean_cum <- rep_curves %>%
  group_by(Watering, SiteID, Days, TidalCat) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE),
    se_prop   = sd(prop, na.rm = TRUE) / sqrt(sum(!is.na(prop))),
    .groups = "drop") %>%
  mutate(lwr = pmax(mean_prop - 1.96 * se_prop, 0),upr = pmin(mean_prop + 1.96 * se_prop, 1))

p_spaghetti_flip <- ggplot() +
  #geom_line(data = rep_curves,aes(x = Days, y = prop, group = IDTag, color = SiteID),
 #   alpha = 0.22,linewidth = 0.6) +
  geom_ribbon(data = mean_cum,aes(x = Days, ymin = lwr, ymax = upr, fill = SiteID),
    alpha = 0.18,color = NA) +
  geom_line(data = mean_cum,aes(x = Days, y = mean_prop, color = SiteID, group = SiteID),
    linewidth = 1.2) +
  facet_wrap(~ Watering, ncol = 2) +
  scale_y_continuous(limits = c(0, 0.3)) +
  labs(x = "Day", y = "Cumulative germination (proportion)",
    title = "Cumulative germination trajectories by watering treatment") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right")

p_spaghetti_flip




## =========================================================
## NEW: Convert alive counts to cumulative total germination
## =========================================================

# Starting from fllong, compute cumulative germination per replicate
fllong_cumgerm <- fllong %>%
  arrange(IDTag, Days) %>%
  group_by(IDTag, Seed.Source, Watering, Seeds) %>%
  mutate(
    alive_count = value,
    cum_germ    = cummax(alive_count),   # running max = total germinated so far
    failures    = Seeds - cum_germ
  ) %>%
  ungroup() %>%
  filter(!is.na(cum_germ), Seeds > 0, cum_germ >= 0, failures >= 0)

## -------------------------
## A) Plot cumulative TOTAL germination curves
## -------------------------
plot_cum_total <- fllong_cumgerm %>%
  mutate(prop = cum_germ / Seeds) %>%
  group_by(Watering, Seed.Source, Days) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    se_prop   = sd(prop, na.rm = TRUE) / sqrt(sum(!is.na(prop))),
    .groups   = "drop"
  ) %>%
  mutate(
    lwr = pmax(mean_prop - 1.96 * se_prop, 0),
    upr = pmin(mean_prop + 1.96 * se_prop, 1)
  )

p_cum_total <- ggplot(plot_cum_total,
                      aes(x = Days, y = mean_prop, color = Watering, group = Watering)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Watering),
              alpha = 0.18, color = NA) +
  facet_wrap(~ Seed.Source) +
  scale_y_continuous(limits = c(0, 0.4)) +
  scale_fill_manual(values = c("High" = "blue", "Low" = "red")) +
  scale_color_manual(values = c("High" = "blue", "Low" = "red")) +
  labs(x = "Day",
       y = "Cumulative total germination (proportion)",
       title = "Cumulative total germination curves") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

p_cum_total

## -------------------------
## B) Final day = TOTAL germinated, not alive remaining
## -------------------------
final_day_total <- fllong_cumgerm %>%
  group_by(IDTag, Seed.Source, Watering) %>%
  slice_max(order_by = Days, n = 1, with_ties = FALSE) %>%
  ungroup()

mod_final_total <- glmmTMB(
  cbind(cum_germ, failures) ~ Seed.Source * Watering,
  data = final_day_total,
  family = betabinomial()
)

summary(mod_final_total)

emm_final_total <- as.data.frame(
  emmeans(mod_final_total, ~ Seed.Source * Watering, type = "response"))

emm_final_total <- emm_final_total %>%
  mutate(
    prob = dplyr::coalesce(
      if ("prob" %in% names(.)) prob else NA_real_,
      if ("response" %in% names(.)) response else NA_real_,
      if ("rate" %in% names(.)) rate else NA_real_
    ),
    lwr = dplyr::coalesce(
      if ("lower.CL" %in% names(.)) lower.CL else NA_real_,
      if ("asymp.LCL" %in% names(.)) asymp.LCL else NA_real_,
      if ("lower.HPD" %in% names(.)) lower.HPD else NA_real_
    ),
    upr = dplyr::coalesce(
      if ("upper.CL" %in% names(.)) upper.CL else NA_real_,
      if ("asymp.UCL" %in% names(.)) asymp.UCL else NA_real_,
      if ("upper.HPD" %in% names(.)) upper.HPD else NA_real_
    )
  )

## Final-day EMMEANS plot using total germination
p_final_emm_total <- ggplot(emm_final_total,
                             aes(x = Seed.Source, y = prob, color = Watering, group = Watering)) +
  geom_point(position = position_dodge(width = 0.35), size = 2) +
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                position = position_dodge(width = 0.35), width = 0.15) +
  geom_line(position = position_dodge(width = 0.35)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  scale_color_manual(values = c("High" = "blue", "Low" = "red")) +
  labs(x = "Seed source",
       y = "Predicted final total germination (proportion)",
       title = "Final-day total germination EMMEANS") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

p_final_emm_total

## Bar version using total germination
p_final_bar_total <- ggplot(emm_final_total,
                            aes(x = Seed.Source, y = prob, fill = Watering)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                position = position_dodge(width = 0.75), width = 0.2) +
  scale_y_continuous(limits = c(0, 0.3)) +
  scale_fill_manual(values = c("High" = "blue", "Low" = "red")) +
  labs(x = "Seed source",
       y = "Predicted final total germination (proportion)",
       title = "Final-day total germination predictions by treatment") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

p_final_bar_total

## -------------------------
## C) Flipped cumulative total germination plot
## -------------------------
mean_cum_total <- fllong_cumgerm %>%
  mutate(prop = cum_germ / Seeds) %>%
  group_by(Watering, Seed.Source, Days) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    se_prop   = sd(prop, na.rm = TRUE) / sqrt(sum(!is.na(prop))),
    .groups = "drop"
  ) %>%
  mutate(
    lwr = pmax(mean_prop - 1.96 * se_prop, 0),
    upr = pmin(mean_prop + 1.96 * se_prop, 1)
  )

p_spaghetti_flip_total <- ggplot() +
  geom_ribbon(data = mean_cum_total,
              aes(x = Days, ymin = lwr, ymax = upr, fill = Seed.Source),
              alpha = 0.18, color = NA) +
  geom_line(data = mean_cum_total,
            aes(x = Days, y = mean_prop, color = Seed.Source, group = Seed.Source),
            linewidth = 1.2) +
  facet_wrap(~ Watering, ncol = 2) +
  scale_y_continuous(limits = c(0, 0.3)) +
  labs(x = "Day",
       y = "Cumulative total germination (proportion)",
       title = "Cumulative total germination by watering treatment") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

p_spaghetti_flip_total
