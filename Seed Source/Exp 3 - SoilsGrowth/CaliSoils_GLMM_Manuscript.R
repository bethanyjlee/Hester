rm(list =ls())

###### GLM Soil Effects on Different Seed Sources and their Germination/Establishment

setwd("~/Hester/Seed Source/Exp 3 - SoilsGrowth")

library(dplyr)
library(tidyverse)



#### Data ####
fulldata <- read.csv("FullGermElkhornSeedSource.csv")

# Reshape from wide (dates) → long
data_long <- fulldata %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Date",
    values_to = "value"
  )

# Clean + format
dataready <- data_long %>%
  mutate(
    Soils = factor(Soil.Type, levels = c("Hester", "Combo", "Potting Soil")),
    Seed.Source = factor(Seed.Source),
    SourceID = factor(SiteID),
    Tidal = factor(TidalCat),
    
    # FIX DATE FORMAT
    Date = gsub("^X", "", Date),        # remove leading X
    Date = gsub("\\.", "/", Date),      # convert . to /
    Date = as.Date(Date, format = "%m/%d/%Y"),
    
    successes = value,
    failures  = Seeds - value
  ) %>%
  filter(Seeds > 0)


data_summary <- dataready %>%
  group_by(TagID) %>%
  filter(Date == max(Date)) %>%
  ungroup() %>%
  mutate(germ_prop = successes / Seeds) %>%
 group_by(SiteID) %>%
  summarise(
    n         = n(),
    mean_germ = mean(germ_prop, na.rm = TRUE),
    sd_germ   = sd(germ_prop, na.rm = TRUE),
    se_germ   = sd_germ / sqrt(n)
  )

print(data_summary)

data_summary_tidal <- dataready %>%
  group_by(TagID) %>%                    # FIX HERE
  filter(Date == max(Date)) %>%
  ungroup() %>%
  mutate(germ_prop = successes / Seeds) %>%
  group_by(Tidal) %>%
  summarise(
    n         = n(),
    mean_germ = mean(germ_prop, na.rm = TRUE),
    sd_germ   = sd(germ_prop, na.rm = TRUE),
    se_germ   = sd_germ / sqrt(n))

print(data_summary_tidal)

#### Model Discovery - By Site ####



# Model with interactions
m_full <- glm(
  cbind(successes, failures) ~ Soils * Seed.Source,
  data = dataready,
  family = binomial)
summary(m_full)

# Model without interactions

m_no_interaction <- glm(
  cbind(successes, failures) ~ Soils + Seed.Source,
  data = dataready,
  family = binomial)

anova(m_full, m_no_interaction, test = "Chisq") # AIC is lower for with interactions
# use m_full but lets check family use

#quasi 
m_quasi <- glm(
  cbind(successes, failures) ~ Soils * Seed.Source,
  data = dataready,
  family = quasibinomial
)
summary(m_quasi)
car::Anova(m_quasi, type = "II")


#betabinomial
library(glmmTMB)
m_bb <- glmmTMB(
  cbind(successes, failures) ~ Soils * SiteID,
  data = dataready,
  family = betabinomial())
summary(m_bb)

#### comparing all AICs, the beta binomial is the best model






#### Model Discovery - Tidal Cat ####


# Model with interactions
m_full_tid <- glm(
  cbind(successes, failures) ~ Soils * Tidal,
  data = dataready,
  family = binomial)
summary(m_full_tid)

# Model without interactions

m_no_interaction_tid <- glm(
  cbind(successes, failures) ~ Soils + Tidal,
  data = dataready,
  family = binomial)

anova(m_full_tid, m_no_interaction_tid, test = "Chisq") # AIC is lower for with interactions
# use m_full but lets check family use

#quasi 
m_quasi_tid <- glm(
  cbind(successes, failures) ~ Soils * Tidal,
  data = dataready,
  family = quasibinomial)
summary(m_quasi_tid)
car::Anova(m_quasi, type = "II")


#betabinomial
library(glmmTMB)
m_bb_tid <- glmmTMB(
  cbind(successes, failures) ~ Soils * Tidal,
  data = dataready,
  family = betabinomial())
summary(m_bb_tid)

fixed_tidal<- glmmTMB(cbind(successes, failures) ~ Soils * Tidal + (1|SiteID),
  family = betabinomial, data = dataready)
summary(fixed_tidal)




# Marginal means for each Soil × Source combo (response scale)
library(emmeans)

emm <- emmeans(fixed_tidal, ~ Soils | Tidal, type = "response")
summary(emm)

# All pairwise combo × combo (Tukey)
pairs(emm, adjust = "tukey")  # report this table

# Within-soil source comparisons
pairs(emm, by = "Soils", adjust = "tukey")

# Within-TidalCat soil comparisons
pairs(emm, by = "Tidal", adjust = "tukey")

# compact letters for figures / tables
library(multcomp)  # for CLD letters
cld(emm, adjust = "tukey", Letters = letters, type = "response")


##### Post hoc analyses ####

library(dplyr)
#### OVERALL MODEL TESTS ####

library(car)

# Type III Wald chi-square tests
car::Anova(fixed_tidal, type = "III")

# pseudo-R2
library(performance)
r2(fixed_tidal)

# model coefficients with CIs
confint(fixed_tidal)

# likelihood ratio tests for main effects + interaction
fixed_no_int <- glmmTMB(
  cbind(successes, failures) ~ Soils + Tidal + (1|SiteID),
  family = betabinomial,
  data = dataready)

anova(fixed_no_int, fixed_tidal)

#### EMMEANS + ODDS RATIOS ####

library(emmeans)

emm_tidal <- emmeans(
  fixed_tidal,
  ~ Soils | Tidal,
  type = "response")

summary(emm_tidal)

# pairwise soil comparisons within tidal category
soil_pairs <- pairs(
  emm_tidal,
  by = "Tidal",
  adjust = "tukey")

summary(soil_pairs, infer = TRUE)

# odds ratios instead of log odds
soil_OR <- contrast(
  emm_tidal,
  method = "pairwise",
  by = "Tidal")

summary(soil_OR, type = "response")

#### SITE-LEVEL EFFECTS ####

fixed_site <- glmmTMB(
  cbind(successes, failures) ~ Soils * SiteID,
  family = betabinomial,
  data = dataready)

car::Anova(fixed_site, type = "III")

emm_site <- emmeans(
  fixed_site,
  ~ SiteID,
  type = "response")

pairs(emm_site, adjust = "tukey")





##### plot graphs prep#####

#### ESTIMATED GERMINATION PERCENTAGES ####

emm_percent <- summary(emm_tidal) %>%
  mutate(percent = prob * 100)

emm_percentlibrary(stringr)
library(ggplot2)
library(scales)
library(multcomp)  # for cld

# Tidy emmeans + letters (response scale)
# 1) Base tables for post hoc plots
emm_df <- as.data.frame(summary(emm)) %>%
  dplyr::mutate(
    Soils = factor(Soils, levels = c("Hester","Combo","Potting Soil")),
    pct   = prob * 100)

letters_df <- multcomp::cld(emm, adjust = "sidak", Letters = letters, type = "response") %>%
  as.data.frame() %>%
  dplyr::transmute(Soils, Tidal, letters = stringr::str_trim(.group))

# 2) Choose an ordering for Seed.Source
order_by_overall <- emm_df %>%
  dplyr::group_by(Tidal) %>%
  dplyr::summarise(m = mean(prob, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(m) %>%
  dplyr::pull(Tidal)

# (Alt) Order by Hester-only performance:
# order_by_hester <- emm_df %>%
#   dplyr::filter(Soils == "Hester") %>%
#   dplyr::arrange(prob) %>% dplyr::pull(Seed.Source)

# 3) Build plot_df safely, then plot
plot_df <- emm_df %>%
  dplyr::left_join(letters_df, by = c("Soils","Tidal")) %>%
  dplyr::mutate(Tidal = factor(Tidal, levels = order_by_overall))

## got to add days after sowing first
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)


sowing_date <- as.Date("2024-03-23")

# Rebuild `long` with robust date parsing from the pivoted header
long <- fulldata %>%
  pivot_longer(starts_with("X"), names_to = "Date", values_to = "value") %>%
  mutate(
    # strip the leading X and normalize separators
    Date_str = gsub("^X", "", Date),
    # parse against multiple possible formats (3/23/24, 03.23.24, 3-23-2024, etc.)
    Date = as.Date(parse_date_time(Date_str,
                                   orders = c("mdy","m.d.y","m/d/y","m-d-y",
                                              "mdY","m.d.Y","m/d/Y","m-d-Y"))),
    value = as.numeric(value)
  ) %>%
  filter(!is.na(Date))

# For plotting only:
dataready <- long %>%
  mutate(
    Soils       = factor(Soil.Type, levels = c("Hester","Combo","Potting Soil")),
    Seed.Source = factor(SiteID),
    Seeds       = as.numeric(Seeds),
    DaysAfterSowing = as.integer(difftime(Date, sowing_date, units = "days")),  # 0 on sow date
    prop = pmin(pmax(value / Seeds, 0), 1),
    Tidal       = factor(TidalCat))






# Observed per-date proportions
raw_ts <- dataready %>%
  mutate(
    DaysAfterSowing = as.integer(Date - sowing_date),
    prop = pmin(pmax(value / Seeds, 0), 1)
  ) %>%
  group_by(DaysAfterSowing, Soils, Seed.Source, TidalCat, SiteID) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    se        = sd(prop,  na.rm = TRUE) / sqrt(n()),
    n         = n(),
    .groups   = "drop"
  )

group_colors<- c(Natural = "#1E88E5", Diked = '#D81B60', Restored = '#FFC107')
neworder<-c("Natural","Diked","Restored")
raw_ts<-arrange(transform(raw_ts, TidalCat = factor(TidalCat, levels = neworder)), TidalCat)


# Plot: lines + points; error bars only every 7 days
gg_raw_ts <- ggplot(raw_ts, aes(x = DaysAfterSowing, y = mean_prop, 
                                fill = TidalCat, color = TidalCat)) +
  geom_smooth(linewidth = 0.7) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_errorbar(
    data = subset(raw_ts, DaysAfterSowing %% 7 == 0),
    aes(ymin = mean_prop - se, ymax = mean_prop + se),
    width = 0.4, alpha = 0.6, color = group_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1), name = "Observed germination (%)") +
  labs(x = "Days after sowing", color = "Seed source") +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors, guide = 'none') +
  facet_grid(~ Soils) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

gg_raw_ts

soil_colors<- c(Hester = "#999333", Combo = '#aa4499', "Potting Soil" = '#661100')

#### MS FIgure 5 ####

tidal_labels <- c(
  Natural = "Natural Provenances",
  Diked = "Diked Provenances",
  Restored = "Restored Provenances"
)

gg_raw_soil <- ggplot(raw_ts, aes(x = DaysAfterSowing, y = mean_prop , 
                                fill = Soils, color = Soils)) +
  geom_smooth(linewidth = 0.7) +
 # geom_point(size = 1.5, alpha = 0.8) +
  geom_errorbar(
    data = subset(raw_ts, DaysAfterSowing %% 7 == 0),
    aes(ymin = mean_prop - se, ymax = mean_prop + se),
    width = 0.4, alpha = 0.6, color = group_colors) +
 # scale_y_continuous(labels = percent_format(accuracy = 1), name = "Observed germination (%)") +
  labs(x = "Days after sowing", color = "Soil Type") +
  scale_color_manual(values = soil_colors) +
  scale_fill_manual(values = soil_colors, guide = 'none') +
  facet_grid(~ TidalCat, labeller = labeller(TidalCat = tidal_labels)) +
  theme_bw(base_size = 12) +
  labs(y = "Proportion Germinated") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

gg_raw_soil
ggsave(filename = "Figure4.tif", dpi = 300, path = "Figures")


# Helper: last non-NA# Helper: last non-NATidalCat
last_non_na <- function(x) { y <- x[!is.na(x)]; if (length(y)) y[length(y)] else NA_real_ }

# Attach a row id to track replicates through long pivot
long_id <- fulldata %>%
  mutate(.rowid = dplyr::row_number()) %>%
  tidyr::pivot_longer(cols = starts_with("X"), names_to = "Date") %>%
  mutate(
    Date = as.Date(gsub("X", "", Date), format = "%m.%d.%y"),
    value = as.numeric(value),
    Soils = factor(Soil.Type, levels = c("Hester","Combo","Potting Soil")),
    Seed.Source = factor(SiteID),
    TidalCat = factor(TidalCat),
    
  )

last_non_na <- function(x) {
  x <- x[!is.na(x)]
  if(length(x) == 0) return(NA)
  tail(x, 1)
}

# Final (cumulative) count per replicate
final_per_rep <- long_id %>%
  arrange(.rowid, Date) %>%
  group_by(.rowid, Soils, SiteID, Seeds, TidalCat) %>%
  summarise(final = last_non_na(value), .groups = "drop") %>%
  mutate(final = pmax(pmin(final, Seeds), 0))

# Aggregate to Soil × Source totals and compute Wilson CI
agg_base <- final_per_rep %>%
  group_by(Soils, TidalCat, SiteID) %>%
  summarise(success = sum(final, na.rm = TRUE),
            trials  = sum(Seeds, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(
    p  = success / trials,
    z  = qnorm(0.975),
    den = 1 + (z^2)/trials,
    ctr = p + (z^2)/(2*trials),
    adj = z * sqrt((p*(1 - p)/trials) + (z^2)/(4*trials^2)),
    lcl = pmax(0, (ctr - adj)/den),
    ucl = pmin(1, (ctr + adj)/den)
  )


# Plot as bars (faceted by soil)
gg_base_final <- ggplot(
  agg_base,
  aes(x = SiteID, y = p, fill = Soils)
) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.7,
           color = "grey30") +
  geom_errorbar(
    aes(ymin = lcl, ymax = ucl),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    name = "Final observed germination (%)"
  ) +
  scale_fill_manual(values = soil_colors) +
  labs(
    x = "Seed source",
    fill = "Soil",
    title = "Final-day observed germination (Wilson 95% CI)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

gg_base_final



# Helper: last non-NA
last_non_na <- function(x) { y <- x[!is.na(x)]; if (length(y)) y[length(y)] else NA_real_ }

# Attach a row id to track replicates through long pivot
long_id <- fulldata %>%
  mutate(.rowid = dplyr::row_number()) %>%
  tidyr::pivot_longer(cols = starts_with("X"), names_to = "Date") %>%
  mutate(
    Date = as.Date(gsub("X", "", Date), format = "%m.%d.%y"),
    value = as.numeric(value),
    Soils = factor(Soil.Type, levels = c("Hester","Combo","Potting Soil")),
    Seed.Source = factor(SiteID),
    TidalCat = factor(TidalCat)
  )

# Final (cumulative) count per replicate
final_per_rep <- long_id %>%
  arrange(.rowid, Date) %>%
  group_by(.rowid, Soils, TidalCat, Seeds) %>%
  summarise(final = last_non_na(value), .groups = "drop") %>%
  mutate(final = pmax(pmin(final, Seeds), 0))

# Aggregate to Soil × Source totals and compute Wilson CI
agg_base <- final_per_rep %>%
  group_by(Soils, TidalCat) %>%
  summarise(success = sum(final, na.rm = TRUE),
            trials  = sum(Seeds, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(
    p  = success / trials,
    z  = qnorm(0.975),
    den = 1 + (z^2)/trials,
    ctr = p + (z^2)/(2*trials),
    adj = z * sqrt((p*(1 - p)/trials) + (z^2)/(4*trials^2)),
    lcl = pmax(0, (ctr - adj)/den),
    ucl = pmin(1, (ctr + adj)/den)
  )

# Optional: order sources by overall observed mean
order_sources_obs <- agg_base %>%
  group_by(TidalCat) %>%
  summarise(m = mean(p), .groups = "drop") %>%
  arrange(m) %>% pull(TidalCat)

agg_base <- agg_base %>%
  mutate(Seed.Source = factor(Seed.Source, levels = order_sources_obs))

# Plot as bars (faceted by soil)
gg_base_final <- ggplot(agg_base, aes(x = TidalCat, y = p, fill = Soils)) +
  geom_col(width = 0.7, color = "grey30") +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), name = "Final observed germination (%)") +
  labs(x = "Seed source", fill = "Soil",
       title = "Final-day observed germination (Wilson 95% CI)") +
  facet_wrap(~ Soils, nrow = 1) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())
gg_base_final


#### FINAL OBSERVED GERMINATION BAR GRAPH — BY INDIVIDUAL SITE ####

# Helper: last non-NA value
last_non_na <- function(x) {
  y <- x[!is.na(x)]
  if (length(y)) y[length(y)] else NA_real_
  }

# Add row IDs and pivot long
long_id <- fulldata %>%
  mutate(.rowid = dplyr::row_number()) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Date",
    values_to = "value") %>%
  mutate(
    Date = gsub("^X", "", Date),
    Date = gsub("\\.", "/", Date),
    Date = as.Date(Date, format = "%m/%d/%Y"),
    value = as.numeric(value),
    Soils = factor(
      Soil.Type,
      levels = c("Hester", "Combo", "Potting Soil")),
    Seed.Source = factor(SiteID),
    TidalCat = factor(TidalCat))

# Get FINAL cumulative germination for each replicate
final_per_rep <- long_id %>%
  arrange(.rowid, Date) %>%
  group_by(.rowid, Soils, Seed.Source, Seeds, TidalCat) %>%
  summarise(
    final = last_non_na(value),
    .groups = "drop") %>%
  mutate(final = pmax(pmin(final, Seeds), 0))

# Aggregate by Soil × Individual Site
agg_site <- final_per_rep %>%
  group_by(Soils, Seed.Source, TidalCat) %>%
  summarise(
    success = sum(final, na.rm = TRUE),
    trials  = sum(Seeds, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    p  = success / trials,
    z   = qnorm(0.975),
    den = 1 + (z^2)/trials,
    ctr = p + (z^2)/(2*trials),
    adj = z * sqrt((p * (1 - p) / trials) +
                     (z^2)/(4 * trials^2)),
    lcl = pmax(0, (ctr - adj)/den),
    ucl = pmin(1, (ctr + adj)/den))

# Order sites by overall germination


#### supplemental fig 5 ####
site_order <- c(
  "Muted1",
  "Muted2",
  "Muted3",
  "Muted4",
  "Muted5",
  "Restored1",
  "Tidal1",
  "Tidal2",
  "Tidal3",
  "Tidal4",
  "Tidal5",
  "Tidal6")

agg_site <- agg_site %>%
  mutate(
    Seed.Source = factor(
      Seed.Source,
      levels = site_order))

# Soil colors
soil_colors <- c(
  Hester = "#F8766D",
  Combo = "#00BA38",
  "Potting Soil" = "#619CFF")

# FINAL PLOT
SupplementalFigure5 <- ggplot(
  agg_site,
  aes(x = Seed.Source, y = p, fill = Soils)) +
  geom_col(width = 0.7,color = "grey30") +
  geom_errorbar(aes(ymin = lcl, ymax = ucl),
    width = 0.2) +
  facet_wrap(~ Soils, nrow = 1) +
  scale_fill_manual(values = soil_colors) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    name = "Final observed germination (%)") +
  labs(
    x = "Seed source",
    fill = "Soil",
    title = "Final-day observed germination (Wilson 95% CI)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "SupplementalFigure5.tif", dpi = 300, path = "Figures")
