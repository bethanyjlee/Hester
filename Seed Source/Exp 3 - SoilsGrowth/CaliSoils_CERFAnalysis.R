###### GLM Soil Effects on Different Seed Sources and their Germination/Establishment

#### set wd
setwd("~/R data/Chp2SeedSource/Exp 3 - SoilsGrowth")


#### Libraries ####
library(tidyverse)
library(lme4)
library(performance)   # for model checks
library(DHARMa)        # for residual diagnostics
library(car)           # for Type II/III Anova

#### Data ####
fulldata <- read.csv("FullGermElkhornSeedSource.csv")

# Reshape & clean
dataready <- fulldata %>%
  mutate(
    Soils = factor(Soil.Type, levels = c("Hester", "Combo", "Potting Soil")),
    Seed.Source = factor(Seed.Source),
    SiteId = factor(SiteID),
    successes = value,
    failures = Seeds - value
  ) %>%
  filter(Seeds > 0)

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
  cbind(successes, failures) ~ Soils * Seed.Source,
  data = dataready,
  family = betabinomial()
)
summary(m_bb)

#### comparing all AICs, the beta binomial is the best model

library(emmeans)

# Marginal means for each Soil × Source combo (response scale)
emm <- emmeans(m_bb, ~ Soils * Seed.Source, type = "response")
summary(emm)         # predicted germination %, SE, CI

# All pairwise combo × combo (Tukey)
pairs(emm, adjust = "tukey")  # report this table

# Within-soil source comparisons
pairs(emm, by = "Soils", adjust = "tukey")

# Within-source soil comparisons
pairs(emm, by = "Seed.Source", adjust = "tukey")

# compact letters for figures / tables
library(multcomp)  # for CLD letters
cld(emm, adjust = "tukey", Letters = letters, type = "response")

##### plot graphs #####
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(multcomp)  # for cld

# Tidy emmeans + letters (response scale)
# 1) Base tables for post hoc plots
emm_df <- as.data.frame(summary(emm)) %>%
  dplyr::mutate(
    Soils = factor(Soils, levels = c("Hester","Combo","Potting Soil")),
    pct   = prob * 100
  )

letters_df <- multcomp::cld(emm, adjust = "sidak", Letters = letters, type = "response") %>%
  as.data.frame() %>%
  dplyr::transmute(Soils, Seed.Source, letters = stringr::str_trim(.group))

# 2) Choose an ordering for Seed.Source
order_by_overall <- emm_df %>%
  dplyr::group_by(Seed.Source) %>%
  dplyr::summarise(m = mean(prob, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(m) %>%
  dplyr::pull(Seed.Source)

# (Alt) Order by Hester-only performance:
# order_by_hester <- emm_df %>%
#   dplyr::filter(Soils == "Hester") %>%
#   dplyr::arrange(prob) %>% dplyr::pull(Seed.Source)

# 3) Build plot_df safely, then plot
plot_df <- emm_df %>%
  dplyr::left_join(letters_df, by = c("Soils","Seed.Source")) %>%
  dplyr::mutate(Seed.Source = factor(Seed.Source, levels = order_by_overall))

#### plots ####
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
    Seed.Source = factor(Seed.Source),
    Seeds       = as.numeric(Seeds),
    DaysAfterSowing = as.integer(difftime(Date, sowing_date, units = "days")),  # 0 on sow date
    prop = pmin(pmax(value / Seeds, 0), 1)
  )






# Observed per-date proportions
raw_ts <- dataready %>%
  mutate(
    DaysAfterSowing = as.integer(Date - sowing_date),
    prop = pmin(pmax(value / Seeds, 0), 1)
  ) %>%
  group_by(DaysAfterSowing, Soils, Seed.Source) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    se        = sd(prop,  na.rm = TRUE) / sqrt(n()),
    n         = n(),
    .groups   = "drop"
  )

# Plot: lines + points; error bars only every 7 days
gg_raw_ts <- ggplot(raw_ts, aes(x = DaysAfterSowing, y = mean_prop, color = Seed.Source)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_errorbar(
    data = subset(raw_ts, DaysAfterSowing %% 7 == 0),
    aes(ymin = mean_prop - se, ymax = mean_prop + se),
    width = 0.4, alpha = 0.6
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), name = "Observed germination (%)") +
  labs(x = "Days after sowing", color = "Seed source") +
  facet_grid(~ Soils) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

gg_raw_ts





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
    Seed.Source = factor(Seed.Source)
  )

# Final (cumulative) count per replicate
final_per_rep <- long_id %>%
  arrange(.rowid, Date) %>%
  group_by(.rowid, Soils, Seed.Source, Seeds) %>%
  summarise(final = last_non_na(value), .groups = "drop") %>%
  mutate(final = pmax(pmin(final, Seeds), 0))

# Aggregate to Soil × Source totals and compute Wilson CI
agg_base <- final_per_rep %>%
  group_by(Soils, Seed.Source) %>%
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
  group_by(Seed.Source) %>%
  summarise(m = mean(p), .groups = "drop") %>%
  arrange(m) %>% pull(Seed.Source)

agg_base <- agg_base %>%
  mutate(Seed.Source = factor(Seed.Source, levels = order_sources_obs))

# Plot as bars (faceted by soil)
gg_base_final <- ggplot(agg_base, aes(x = Seed.Source, y = p, fill = Soils)) +
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
    Seed.Source = factor(Seed.Source)
  )

# Final (cumulative) count per replicate
final_per_rep <- long_id %>%
  arrange(.rowid, Date) %>%
  group_by(.rowid, Soils, Seed.Source, Seeds) %>%
  summarise(final = last_non_na(value), .groups = "drop") %>%
  mutate(final = pmax(pmin(final, Seeds), 0))

# Aggregate to Soil × Source totals and compute Wilson CI
agg_base <- final_per_rep %>%
  group_by(Soils, Seed.Source) %>%
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
  group_by(Seed.Source) %>%
  summarise(m = mean(p), .groups = "drop") %>%
  arrange(m) %>% pull(Seed.Source)

agg_base <- agg_base %>%
  mutate(Seed.Source = factor(Seed.Source, levels = order_sources_obs))

# Plot as bars (faceted by soil)
gg_base_final <- ggplot(agg_base, aes(x = Seed.Source, y = p, fill = Soils)) +
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




gg_base_points <- ggplot(agg_base, aes(x = Seed.Source, y = p)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.15) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Seed source", y = "Final observed germination (%)",
       title = "Observed final germination by source within soil") +
  facet_wrap(~ Soils, nrow = 1) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg_base_points









## post hoc plots
gg_heat <- ggplot(plot_df, aes(x = Soils, y = Seed.Source, fill = prob)) +
  geom_tile(color = "white", linewidth = 0.3) +
  # Show % and letters (stacked)
  geom_text(aes(label = sprintf("%.1f%%\n(%s)", pct, letters)), size = 3.25, lineheight = 0.95) +
  scale_fill_viridis_c(name = "Predicted\nGermination", labels = percent_format(accuracy = 1)) +
  labs(x = "Soil", y = "Seed source",
       title = "Predicted germination by Soil × Seed Source (beta-binomial emmeans)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

gg_heat
# ggsave("fig_heatmap_soil_source.png", gg_heat, width = 7, height = 6, dpi = 300)






gg_interact <- ggplot(plot_df,
                      aes(x = Seed.Source, y = prob, color = Soils)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                position = position_dodge(width = 0.5), width = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), name = "Predicted germination (%)") +
  labs(x = "Seed source", color = "Soil",
       title = "Interaction of soil and seed source (emmeans, response scale)") +
  facet_wrap(~ Soils, nrow = 1) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

gg_interact
# ggsave("fig_interaction_facets.png", gg_interact, width = 10, height = 4.5, dpi = 300)

gg_bar <- ggplot(plot_df, aes(x = Seed.Source, y = prob, fill = Soils)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72, color = "grey30") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_text(aes(label = letters, y = asymp.UCL + 0.01),
            position = position_dodge(width = 0.8), size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 1), name = "Predicted germination (%)") +
  labs(x = "Seed source", fill = "Soil",
       title = "Soil × seed source emmeans with compact-letter display") +
  facet_wrap(~ Soils, nrow = 1) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

gg_bar
# ggsave("fig_bars_CLD.png", gg_bar, width = 10, height = 4.8, dpi = 300)
