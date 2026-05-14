###### Wesley Experiment - Buried v Surface Seeds

setwd("~/R data/Chp1PWSeeds/Zea and Wesley Exp")

#### import data ####
dat <- read.csv("WesBurSur.csv")

#### load libraries ####
library(ggplot2)
library(dplyr)
library(tidyr)
library(brglm2)
library(emmeans)

#### ---- CLEAN + PREP ---- ####

# Standardize Soil
dat <- dat %>%
  mutate(
    Soil = as.character(Soil),
    Soil = ifelse(Soil %in% c("H1", "H2"), "Hester", Soil),
    Soil = factor(Soil)
  )

# Identify day columns
day_cols <- grep("^D\\d+$", names(dat), value = TRUE)

# Ensure numeric
dat <- dat %>%
  mutate(across(all_of(day_cols), ~ as.numeric(.)))

# Compute Max Germination if missing
if (!"MaxGerm" %in% names(dat)) {
  dat <- dat %>%
    mutate(MaxGerm = do.call(pmax, c(across(all_of(day_cols)), na.rm = TRUE)))
}

# Assign total seeds
if (!"Total" %in% names(dat)) {
  dat <- dat %>% mutate(Total = 50)
}

#### ---- GLM ---- ####

mod_binom <- glm(cbind(MaxGerm, Total - MaxGerm) ~ Soil + Treatment,
                 family = binomial,
                 data = dat)

mod_binom_int <- glm(cbind(MaxGerm, Total - MaxGerm) ~ Soil * Treatment,
                     family = binomial,
                     data = dat)

AIC(mod_binom, mod_binom_int)

# Bias-reduced model
mod_br <- glm(cbind(MaxGerm, Total - MaxGerm) ~ Soil * Treatment,
              family = binomial,
              data = dat,
              method = "brglmFit")

summary(mod_br)

#### ---- PLOT DATA ---- ####

dat <- dat %>%
  mutate(
    germ_pct = (MaxGerm / Total) * 100
  )

# Clean soil labels
dat <- dat %>%
  mutate(
    Soil = case_when(
      Soil %in% c("Hester", "1") ~ "Hester Soil",
      Soil == "PS" ~ "Potting Soil",
      TRUE ~ as.character(Soil)
    )
  )

dat$Soil <- factor(dat$Soil,
                   levels = c("Hester Soil", "Potting Soil"))

#### ---- PLOT ---- ####

ggplot(dat, aes(x = Treatment, y = germ_pct, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  facet_wrap(~ Soil) +
  labs(x = "Treatment", y = "Germination (%)") +
  scale_y_continuous(limits = c(0, 50)) +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

ggsave("LeeFigure5.tif", dpi=300)
