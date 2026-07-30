rm(list =ls())

# Seed Viability Analysis
# Nonlinear Beta-Binomial Model

library(dplyr)
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(performance)

#### Import Data ####


setwd("~/Hester/Seed Source/Exp 1 - Viability")

viab <- read.csv("Viability.csv") %>%
  mutate(
    Site = factor(SiteID),
    TidalCat = factor(TidalCat),
    Fail = Total - Seeds) %>%
  filter(Total > 0)


# Center Elevation

mean_elev <- mean(viab$Elevation, na.rm = TRUE)

viab <- viab %>%
  mutate(
    Elev_c = Elevation - mean_elev)

#### Summary Table ####

data_summary <- viab %>%
  group_by(TidalCat) %>%
  summarise(
    n = n(),
    mean_viab = mean(Seeds / Total),
    sd_viab = sd(Seeds / Total),
    se_viab = sd_viab / sqrt(n))

print(data_summary)


#### Model Selection ####
 

# Linear model
m_linear <- glmmTMB(
  cbind(Seeds, Fail) ~
    TidalCat * Elev_c,
  family = betabinomial(link = "logit"),
  data = viab)

# Quadratic model
m_quad <- glmmTMB(
  cbind(Seeds, Fail) ~
    TidalCat * (Elev_c + I(Elev_c^2)),
  family = betabinomial(link = "logit"),
  data = viab)


# Compare Models


AIC(m_linear, m_quad)

anova(m_linear, m_quad)

summary(m_quad)


# subset models 
restored <- subset(viab, TidalCat == "Restored")

m_rest_linear <- glmmTMB(
  cbind(Seeds, Fail) ~ Elev_c,
  family = betabinomial(link = "logit"),
  data = restored)

m_rest_quad <- glmmTMB(
  cbind(Seeds, Fail) ~ Elev_c + I(Elev_c^2),
  family = betabinomial(link = "logit"),
  data = restored)

anova(m_rest_linear, m_rest_quad)

# Final Model


m_quad2 <- glmmTMB(
  cbind(Seeds, Fail) ~
    TidalCat * (Elev_c + I(Elev_c^2)),
  family = betabinomial(link = "logit"),
  data = viab)

summary(m_quad2)

car::Anova(m_quad2, type = 3)


#### Estimated Marginal Means ####

emm <- emmeans(
  m_quad2,
  ~ TidalCat,
  type = "response")

print(emm)

pairs(emm)



# Calculate Peak Elevation (Restored Marsh)

coefs <- fixef(m_quad2)$cond

b1_restored <-
  coefs["Elev_c"] +
  coefs["TidalCatRestored:Elev_c"]

b2_restored <-
  coefs["I(Elev_c^2)"] +
  coefs["TidalCatRestored:I(Elev_c^2)"]

peak_centered <- -b1_restored / (2 * b2_restored)

peak_elevation <- mean_elev + peak_centered

cat("\nPeak elevation for restored marshes:\n")
print(peak_elevation)


# Prediction Grid


pred_grid <- expand.grid(
  Elevation = seq(
    min(viab$Elevation),
    max(viab$Elevation),
    length.out = 200),
  TidalCat = levels(viab$TidalCat))

pred_grid$Elev_c <-
  pred_grid$Elevation - mean_elev

pred_link <- predict(
  m_quad2,
  newdata = pred_grid,
  type = "link",
  se.fit = TRUE)

pred_grid <- pred_grid %>%
  mutate(
    fit_link = pred_link$fit,
    se_link = pred_link$se.fit,
    fit = plogis(fit_link),
    lwr = plogis(fit_link - 1.96 * se_link),
    upr = plogis(fit_link + 1.96 * se_link))


# Colors

group_colors <- c(
  Natural = "#1E88E5",
  Restored = "#FFC107")

#### Individual Sites Graphs ####

p_indiv <- ggplot() +
  geom_point(
    data = viab,
    aes(
      Elevation,
      Seeds / Total,
      color = Site),
    alpha = 0.4,
    size = 2) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0.02)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)) +
  labs(
    x = "Elevation (m NAVD88)",
    y = "Proportion viable")

p_indiv


ggsave(
  "SupplementalFigure2.tif",
  p_main,
  dpi = 300,
  width = 7,
  height = 5)

#### R² for each marsh type ####

# Natural predictions
natural_dat <- filter(viab, TidalCat == "Natural")

pred_nat <- predict(
  m_quad2,
  newdata = natural_dat,
  type = "response")

r2_nat <- round(
  cor(
    natural_dat$Seeds / natural_dat$Total,
    pred_nat)^2,3)

# Restored predictions
restored_dat <- filter(viab, TidalCat == "Restored")

pred_rest <- predict(
  m_quad2,
  newdata = restored_dat,
  type = "response")

r2_rest <- round(
  cor(
    restored_dat$Seeds / restored_dat$Total,
    pred_rest)^2,3)

r2_labels <- data.frame(
  TidalCat = c("Natural","Restored"),
  label = c(
    paste0("R² = ", r2_nat),
    paste0("R² = ", r2_rest)),
  x = c(
    min(viab$Elevation) + 0.03,
    min(viab$Elevation) + 0.03),
  y = c(0.95,0.87))


#### Figure 1: Combined Plot ####

p_main <- ggplot() +
  geom_point(
    data = viab,
    aes(
      Elevation,
      Seeds / Total,
      color = TidalCat),
    alpha = 0.4,
    size = 2) +
  geom_ribbon(
    data = pred_grid,
    aes(
      x = Elevation,
      ymin = lwr,
      ymax = upr,
      fill = TidalCat),
    alpha = 0.20,
    color = NA) +
  geom_line(
    data = pred_grid,
    aes(
      Elevation,
      fit,
      color = TidalCat),
    linewidth = 1.3) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0.02)) +
  geom_text(
    data = r2_labels,
    aes(
      x = x,
      y = y,
      label = label,
      color = TidalCat),
    hjust = 0,
    size = 5,
    show.legend = FALSE) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)) +
  labs(
    x = "Elevation (m NAVD88)",
    y = "Proportion viable",
    color = "Marsh Type",
    fill = "Marsh Type")

p_main



ggsave(
  "Figure2.tif",
  p_main,
  dpi = 300,
  width = 7,
  height = 5)



#### Diagnostics ####

check_overdispersion(m_quad2)

