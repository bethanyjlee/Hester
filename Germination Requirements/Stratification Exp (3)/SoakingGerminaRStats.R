##### Cleaned Stratification Statistics Code
##### 9/19/2025

## ---- 0) Packages ----
library(GerminaR)
library(dplyr)
library(tidyr)
library(emmeans)
library(broom)
library(performance)   # r2(), check_overdispersion()
library(DHARMa)        # simulation-based diagnostics
library(pscl)          # pR2() McFadden for binomial

## ---- 1) Load NEW-count data & build GerminaR summary ----
# SeedSoaking: daily NEW germination counts in columns D1..Dk + metadata
new <- read.csv("seedSoaking.csv", check.names = FALSE)

# Ensure day columns are numeric and named with prefix "D"
day_cols <- grep("^D\\d+$", names(new), value = TRUE)
stopifnot(length(day_cols) > 0)

new[day_cols] <- lapply(new[day_cols], function(x) suppressWarnings(as.numeric(x)))

# GerminaR summary at the experimental-unit level
gsm <- ger_summary(SeedN = "seeds", evalName = "D", data = new)

# Harmonize factor structure used in models
gsm <- gsm %>%
  mutate(
    Soil.Type     = factor(Soil.Type),
    Soak.Temp     = factor(Soak.Temp),
    Soak.Time     = factor(Soak.Time,     levels = c("3hr","1 day","7 days","Control")),
    Salinity.Soak = factor(Salinity.Soak, levels = c("Freshwater","15ppt","30ppt","Control"))
  )

# Optional: quick table of means ± SE for reporting
sum_by_trt <- gsm %>%
  group_by(Soil.Type, Soak.Temp, Soak.Time, Salinity.Soak) %>%
  summarise(
    n        = dplyr::n(),
    seedsN   = mean(seeds, na.rm = TRUE),
    grs_mean = mean(grs,   na.rm = TRUE),
    grs_se   = sd(grs,     na.rm = TRUE)/sqrt(n),
    mgt_mean = mean(mgt,   na.rm = TRUE),
    mgt_se   = sd(mgt,     na.rm = TRUE)/sqrt(n),
    .groups  = "drop"
  )

## ---- 2) GLM A: Germination proportion (binomial; grs vs seeds-grs) ----
dat_gp <- gsm %>%
  transmute(
    grs, ng = pmax(seeds - grs, 0L),  # non-germinated
    Soil.Type, Soak.Temp, Soak.Time, Salinity.Soak
  )

# Fit binomial GLM
mod_GP <- glm(
  cbind(grs, ng) ~ Soil.Type + Soak.Temp + Soak.Time * Salinity.Soak,
  data = dat_gp,
  family = binomial
)

# Overdispersion check; refit if needed
od <- check_overdispersion(mod_GP)
od_ratio <- od$dispersion_ratio

if (!is.na(od_ratio) && od_ratio > 1.5) {
  mod_GP <- glm(
    cbind(grs, ng) ~ Soil.Type + Soak.Temp + Soak.Time * Salinity.Soak,
    data = dat_gp,
    family = quasibinomial
  )
  gp_family <- "quasibinomial (logit)"
} else {
  gp_family <- "binomial (logit)"
}


summary(mod_GP)

# Binomial diagnostics
sim_gp <- simulateResiduals(mod_GP, n = 1000)
# plot(sim_gp)  # optional; visually inspect if desired

# Fit metrics
gp_summ <- broom::tidy(mod_GP)
gp_glance <- broom::glance(mod_GP)
gp_r2 <- tryCatch(pscl::pR2(mod_GP), error = function(e) NULL)  # McFadden, etc.
gp_perf <- tryCatch(performance::r2(mod_GP), error = function(e) NULL)  # includes Tjur for binomial

## Post-hoc (Tukey) on soaking interaction
emm_gp <- emmeans(mod_GP, ~ Soak.Time * Salinity.Soak, type = "response")
pairs_gp_time_within_sal <- contrast(emm_gp, by = "Salinity.Soak", method = "tukey")
pairs_gp_sal_within_time <- contrast(emm_gp, by = "Soak.Time",     method = "tukey")

# Main-effects EMMS (if you also want them)
emm_gp_soil <- emmeans(mod_GP, ~ Soil.Type,     type = "response")
emm_gp_temp <- emmeans(mod_GP, ~ Soak.Temp,     type = "response")

## ---- 3) GLM B: Mean Germination Time (Gamma with log link) ----
dat_mgt <- gsm %>%
  filter(is.finite(mgt), mgt > 0) %>%   # Gamma requires positive response
  select(mgt, Soil.Type, Soak.Temp, Soak.Time, Salinity.Soak)

mod_MGT <- glm(
  mgt ~ Soil.Type + Soak.Temp + Soak.Time * Salinity.Soak,
  data = dat_mgt,
  family = Gamma(link = "log")
)


summary(mod_MGT)

# Diagnostics
sim_mgt <- simulateResiduals(mod_MGT, n = 1000)
# plot(sim_mgt)  # optional

# Metrics
mgt_summ <- broom::tidy(mod_MGT, exponentiate = TRUE)  # exp(coef) = multiplicative effect on mgt
mgt_glance <- broom::glance(mod_MGT)
mgt_perf <- tryCatch(performance::r2(mod_MGT), error = function(e) NULL)

## Post-hoc (Tukey) on soaking interaction for MGT
emm_mgt <- emmeans(mod_MGT, ~ Soak.Time * Salinity.Soak, type = "response")
pairs_mgt_time_within_sal <- contrast(emm_mgt, by = "Salinity.Soak", method = "tukey")
pairs_mgt_sal_within_time <- contrast(emm_mgt, by = "Soak.Time",     method = "tukey")

# Optional: simple main-effects EMMS
emm_mgt_soil <- emmeans(mod_MGT, ~ Soil.Type, type = "response")
emm_mgt_temp <- emmeans(mod_MGT, ~ Soak.Temp, type = "response")

## ---- 4) (Optional) Soil-stratified models (Potting vs Hester Mimic) ----
fit_by_soil <- function(soil_level) {
  d_gp  <- dat_gp  %>% filter(Soil.Type == soil_level) %>% droplevels()
  d_mgt <- dat_mgt %>% filter(Soil.Type == soil_level) %>% droplevels()
  
  m_gp <- glm(cbind(grs, ng) ~ Soak.Temp + Soak.Time * Salinity.Soak,
              data = d_gp, family = binomial)
  if (check_overdispersion(m_gp)$dispersion_ratio > 1.5) {
    m_gp <- glm(cbind(grs, ng) ~ Soak.Temp + Soak.Time * Salinity.Soak,
                data = d_gp, family = quasibinomial)
  }
  m_mgt <- glm(mgt ~ Soak.Temp + Soak.Time * Salinity.Soak,
               data = d_mgt, family = Gamma(link = "log"))
  list(gp = m_gp, mgt = m_mgt)
}

# Example:
fits_potting <- fit_by_soil("Potting Soil")
fits_hester  <- fit_by_soil("Hester Mimic")

## ---- 5) Compact objects to print/export in your script ----
# Model summaries:
mod_GP; gp_family; gp_glance; gp_r2; gp_perf
mod_MGT; mgt_glance; mgt_perf

# EMMeans tables (use summary() to view with CIs and p-values):
emm_gp; pairs_gp_time_within_sal; pairs_gp_sal_within_time
emm_mgt; pairs_mgt_time_within_sal; pairs_mgt_sal_within_time

# Quick compact coefficient tables (ready to kable/gt if you like):
coef_GP  <- broom::tidy(mod_GP, conf.int = TRUE)
coef_MGT <- broom::tidy(mod_MGT, conf.int = TRUE, exponentiate = TRUE)  # multiplicative effects

# ---- helpers ----
safe_AIC <- function(m) { suppressWarnings(tryCatch(AIC(m), error = function(e) NA_real_)) }
safe_Tjur <- function(m) {
  out <- tryCatch(performance::r2(m), error = function(e) NULL)
  if (is.null(out)) return(NA_real_)
  # performance::r2() returns a named vector or data frame depending on class
  vals <- suppressWarnings(as.data.frame(out))
  if ("R2_Tjur" %in% names(vals)) return(vals$R2_Tjur[1])
  if ("R2" %in% names(vals)) return(vals$R2[1])  # fallback
  NA_real_
}
safe_McFadden <- function(m) {
  out <- tryCatch(pscl::pR2(m), error = function(e) NULL)
  if (is.null(out)) return(NA_real_)
  as.numeric(out[["McFadden"]])
}
safe_R2 <- function(m) {
  out <- tryCatch(performance::r2(m), error = function(e) NULL)
  if (is.null(out)) return(NA_real_)
  # grab first numeric R2 reported
  vals <- suppressWarnings(as.data.frame(out))
  rcols <- names(vals)[sapply(vals, is.numeric)]
  if (length(rcols)) return(vals[[rcols[1]]][1])
  NA_real_
}

# ---- 6A) Germination proportion (grs) model outputs ----
gp_coef <- broom::tidy(mod_GP, conf.int = TRUE) %>%
  mutate(
    odds_ratio = exp(estimate),
    OR_low = exp(conf.low),
    OR_high = exp(conf.high)
  )

gp_fit <- tibble(
  model = "Germination proportion (grs)",
  family = if (exists("gp_family")) gp_family else as.character(family(mod_GP)$family),
  link   = as.character(family(mod_GP)$link),
  AIC    = safe_AIC(mod_GP),             # NA for quasibinomial (expected)
  R2_Tjur = safe_Tjur(mod_GP),
  R2_McFadden = safe_McFadden(mod_GP)
)

# EMMeans on Soak.Time × Salinity.Soak (response scale)
emm_gp_tbl <- as.data.frame(summary(emm_gp, infer = c(TRUE, TRUE))) %>%
  rename(
    mean_resp = response,
    SE = SE,
    df = df,
    lower.CL = lower.CL,
    upper.CL = upper.CL
  )

pairs_gp_time_in_sal <- as.data.frame(summary(pairs_gp_time_within_sal, infer = c(TRUE, TRUE))) %>%
  mutate(scope = "Within salinity: time contrasts")

pairs_gp_sal_in_time <- as.data.frame(summary(pairs_gp_sal_within_time, infer = c(TRUE, TRUE))) %>%
  mutate(scope = "Within time: salinity contrasts")

# ---- 6B) Mean germination time (mgt; Gamma log) outputs ----
mgt_coef <- broom::tidy(mod_MGT, conf.int = TRUE, exponentiate = TRUE) %>%
  rename(
    ratio_mgt = estimate,
    CI_low = conf.low,
    CI_high = conf.high
  )
mgt_fit <- tibble(
  model = "Mean germination time (mgt)",
  family = as.character(family(mod_MGT)$family),
  link   = as.character(family(mod_MGT)$link),
  AIC    = safe_AIC(mod_MGT),
  R2_any = safe_R2(mod_MGT)   # typically Nagelkerke or related for GLMs
)

emm_mgt_tbl <- as.data.frame(summary(emm_mgt, infer = c(TRUE, TRUE))) %>%
  rename(
    mean_mgt = response,
    SE = SE,
    df = df,
    lower.CL = lower.CL,
    upper.CL = upper.CL
  )

pairs_mgt_time_in_sal <- as.data.frame(summary(pairs_mgt_time_within_sal, infer = c(TRUE, TRUE))) %>%
  mutate(scope = "Within salinity: time contrasts")

pairs_mgt_sal_in_time <- as.data.frame(summary(pairs_mgt_sal_within_time, infer = c(TRUE, TRUE))) %>%
  mutate(scope = "Within time: salinity contrasts")

# ---- 6C) (Optional) Soil-stratified fits, if you ran fit_by_soil() above ----
extract_by_soil <- function(fits, soil_label) {
  if (is.null(fits)) return(NULL)
  gp <- fits$gp; mgt <- fits$mgt
  tibble(
    soil = soil_label,
    gp_family = as.character(family(gp)$family),
    gp_link   = as.character(family(gp)$link),
    gp_AIC    = safe_AIC(gp),
    gp_R2_Tjur = safe_Tjur(gp),
    gp_R2_McFadden = safe_McFadden(gp),
    mgt_family = as.character(family(mgt)$family),
    mgt_link   = as.character(family(mgt)$link),
    mgt_AIC    = safe_AIC(mgt),
    mgt_R2_any = safe_R2(mgt)
  )
}
by_soil_summary <- NULL
if (exists("fits_potting")) by_soil_summary <- bind_rows(by_soil_summary, extract_by_soil(fits_potting, "Potting Soil"))
if (exists("fits_hester"))  by_soil_summary <- bind_rows(by_soil_summary, extract_by_soil(fits_hester,  "Hester Mimic"))

# ---- 6D) Console previews (for quick copy into Results text) ----
cat("\n=== MODEL FITS ===\n"); print(gp_fit); print(mgt_fit)
cat("\n=== SIGNIFICANCE (grs) — coefficients ===\n"); print(gp_coef)
cat("\n=== SIGNIFICANCE (mgt) — coefficients (ratios, log link) ===\n"); print(mgt_coef)
cat("\n=== EMMeans (grs; prop germinated) by Soak.Time × Salinity.Soak ===\n"); print(head(emm_gp_tbl, 12))
cat("\n=== EMMeans (mgt) by Soak.Time × Salinity.Soak ===\n"); print(head(emm_mgt_tbl, 12))
cat("\n=== Tukey contrasts (grs) ===\n"); print(head(pairs_gp_time_in_sal, 12)); print(head(pairs_gp_sal_in_time, 12))
cat("\n=== Tukey contrasts (mgt) ===\n"); print(head(pairs_mgt_time_in_sal, 12)); print(head(pairs_mgt_sal_in_time, 12))
if (!is.null(by_soil_summary)) { cat("\n=== Soil-stratified model summaries ===\n"); print(by_soil_summary) }



