## CUMULATIVE TOTAL GERMINATION ANALYSIS

setwd("~/Hester/Seed Source/Exp 2 - Moisture")

library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(cowplot)

#### 1. Read data and combine seed sources ####
fldata <- read.csv("MoistureStress.csv")

fldata <- fldata %>%
  mutate(
    Seed.Source = as.character(Seed.Source),
    Seed.Source = ifelse(Seed.Source %in% c("P1", "P2"), "Hester", Seed.Source),
    Seed.Source = factor(Seed.Source),
    SiteID = factor(SiteID),
    Watering = factor(Watering, levels = c("High", "Low")),
    TidalCat = factor(TidalCat)
  )

#### 2. Reshape to long format ####
fllong <- fldata %>%
  pivot_longer(
    cols = starts_with("D"),
    names_to = "Days",
    values_to = "value"
  ) %>%
  mutate(
    Days = as.numeric(gsub("D", "", Days))
  )

#### 3. Convert alive counts to cumulative germination ####
cum_dat <- fllong %>%
  arrange(IDTag, Days) %>%
  group_by(IDTag, SiteID, Watering, Seeds, Seed.Source, TidalCat) %>%
  mutate(
    alive_count = value,
    cum_germ    = cummax(alive_count),
    failures    = Seeds - cum_germ,
    prop        = cum_germ / Seeds
  ) %>%
  ungroup() %>%
  filter(!is.na(cum_germ), Seeds > 0, cum_germ >= 0, failures >= 0)

#### 4. Option 1: Proportion data ####
plot_dat_prop <- cum_dat %>%
  group_by(SiteID, Seed.Source, Watering, Days, TidalCat) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    se_prop   = sd(prop, na.rm = TRUE) / sqrt(sum(!is.na(prop))),
    .groups   = "drop"
  ) %>%
  mutate(
    lwr = pmax(mean_prop - 1.96 * se_prop, 0),
    upr = pmin(mean_prop + 1.96 * se_prop, 1)
  )

#### 4. Option 2: Count data ####
plot_dat_count <- cum_dat %>%
  group_by(SiteID, Seed.Source, Watering, Days, TidalCat) %>%
  summarise(
    mean_germ = mean(cum_germ, na.rm = TRUE),
    se_germ   = sd(cum_germ, na.rm = TRUE) / sqrt(sum(!is.na(cum_germ))),
    .groups   = "drop"
  ) %>%
  mutate(
    lwr = pmax(mean_germ - 1.96 * se_germ, 0),
    upr = mean_germ + 1.96 * se_germ
  )

#### 5. Set site order by tidal category ####
tidal_sites  <- sort(unique(plot_dat_count$SiteID[plot_dat_count$TidalCat == "Natural"]))
muted_sites  <- sort(unique(plot_dat_count$SiteID[plot_dat_count$TidalCat == "Diked"]))
rest_sites   <- sort(unique(plot_dat_count$SiteID[plot_dat_count$TidalCat == "Restored"]))

panel_order <- c(muted_sites, tidal_sites, rest_sites)

plot_dat_count$SiteID <- factor(plot_dat_count$SiteID, levels = panel_order)
plot_dat_prop$SiteID  <- factor(plot_dat_prop$SiteID,  levels = panel_order)

#### 6.  Count plot ####
watering_colors <- c(Low = 'orange', High = 'skyblue')

final_count_plot <- ggplot(plot_dat_count,
                           aes(x = Days, y = mean_germ,
                               color = Watering, group = Watering)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Watering),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ SiteID, nrow = 2) +
  scale_color_manual(values = watering_colors) +
  scale_fill_manual(values = watering_colors) +
  labs(x = "Day", y = "Cumulative Germination (counts)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

final_count_plot

#### 8. Proportion plot ####
final_cum_plot <- ggplot(plot_dat_prop,
                         aes(x = Days, y = mean_prop,
                             color = Watering, group = Watering)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Watering),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ SiteID, nrow = 2, drop = FALSE) +
  scale_color_manual(values = watering_colors) +
  scale_fill_manual(values = watering_colors) +
  scale_y_continuous(limits = c(0, 0.3)) +
  labs(x = "Day", y = "Proportion Germinated") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

final_cum_plot

#### 9. Stats at final timepoint####

final_day <- cum_dat %>%
  group_by(IDTag, Seed.Source, SiteID, Watering, TidalCat) %>%
  slice_max(order_by = Days, n = 1, with_ties = FALSE) %>%
  ungroup()

stat_dat <- final_day %>%
  group_by(SiteID) %>%
  summarise(
    p_value = {
      m <- glm(cbind(cum_germ, failures) ~ Watering,
               family = binomial, data = cur_data())
      summary(m)$coefficients["WateringLow", "Pr(>|z|)"]
    },
    .groups = "drop"
  ) %>%
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      #  TRUE            ~ "ns"
    )
  )

# get y positions
y_pos <- plot_dat_prop %>%
  group_by(SiteID) %>%
  summarise(y = max(upr, na.rm = TRUE), .groups = "drop")

stat_dat <- left_join(stat_dat, y_pos, by = "SiteID")

#### 10. Add stats to proportion plot ####

final_cum_plot +
  geom_text(data = stat_dat,
            aes(x = Inf, y = y, label = sig),
            inherit.aes = FALSE,
            hjust = 1.2,
            size = 5,
            fontface = "bold")

final_cum_plot +geom_text(data = stat_dat,aes(x = Inf, y = y, label = sig),inherit.aes = FALSE,hjust = 1.2,size = 5,fontface = "bold")

ggsave(filename = "SupplementalFig3.tif", dpi = 300, path = "Figures")                                                                                                                                                                                                                                                   
