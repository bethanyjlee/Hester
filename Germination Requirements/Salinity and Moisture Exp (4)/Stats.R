####### Inhibition Experiment Stats
#### 9/19/2025

#### 0) Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(GerminaR)
  library(car)
  library(emmeans)
  library(DHARMa)
  library(performance)
  library(glmmTMB)
  library(broom)
})

#### 1) Load & prepare raw daily data ------------------------------------------
# Point this to your file with daily counts (NEW or cumulative): D1, D2, ... DN
file_germ <- "Mega2New.csv"
raw <- read.csv(file_germ, check.names = FALSE)

# Identify day columns safely
day_cols <- grep("^D\\d+$", names(raw), value = TRUE)
if (length(day_cols) == 0) stop("No day columns like D1..D* found in the file.")

to_num <- function(x) suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(x)), fixed = TRUE)))

dat <- raw %>%
  mutate(
    across(all_of(day_cols), ~ pmax(to_num(.), 0)),
    Seeds     = as.integer(Seeds),
    Soil      = factor(Soil),
    Salinity  = factor(Salinity),
    Moisture  = factor(Moisture)
  ) %>%
  drop_na(Seeds, Soil, Salinity, Moisture) %>%
  filter(Seeds > 0)

# Final germinated per unit (robust to daily-new vs cumulative inputs)
daily_sum <- dat %>% select(all_of(day_cols)) %>% as.matrix() %>% rowSums(na.rm = TRUE)
last_val  <- dat %>% select(all_of(day_cols)) %>% as.matrix() %>% apply(1, function(r){
  x <- r[!is.na(r)]; if (length(x)==0) return(0); tail(x,1)
})
final_candidate <- pmax(daily_sum, last_val, na.rm = TRUE)
final_germ <- pmin(pmax(round(final_candidate), 0), dat$Seeds)

final <- dat %>%
  mutate(
    Success = final_germ,
    Failure = pmax(Seeds - final_germ, 0)
  ) %>%
  select(Soil, Salinity, Moisture, Seeds, Success, Failure, all_of(day_cols))

# Optional: set factor orders for manuscript clarity
final <- final %>%
  mutate(
    Moisture = factor(Moisture, levels = c("Low","Medium","High")),
    Salinity = fct_relevel(Salinity, sort(levels(Salinity))),
    Soil     = fct_relevel(Soil, sort(levels(Soil)))
  )

#### 2) GERMINATION PERCENTAGE MODEL -------------------------------------------
cat("\n================ GERMINATION PERCENTAGE MODEL ================\n")

# Binomial (exploratory; confirm overdispersion)
m_bin_full3 <- glm(cbind(Success, Failure) ~ Soil * Salinity * Moisture,
                   data = final, family = binomial)

cat("\n[Binomial] Overdispersion check:\n")
print(performance::check_overdispersion(m_bin_full3))

# Beta-binomial (inference model)
m_bb_main  <- glmmTMB(cbind(Success, Failure) ~ Soil + Salinity + Moisture,
                      data = final, family = betabinomial(link = "logit"))
m_bb_2way  <- glmmTMB(cbind(Success, Failure) ~ (Soil + Salinity + Moisture)^2,
                      data = final, family = betabinomial(link = "logit"))
m_bb_full3 <- glmmTMB(cbind(Success, Failure) ~ Soil * Salinity * Moisture,
                      data = final, family = betabinomial(link = "logit"))

cat("\n[Beta-binomial] Model comparisons (LR tests):\n")
cat("  2-way vs main effects:\n")
print(anova(m_bb_main,  m_bb_2way))
cat("  3-way vs 2-way:\n")
print(anova(m_bb_2way,  m_bb_full3))

cat("\nAIC progression (lower is better):\n")
fit_tab <- tibble(
  model  = c("binomial_full3","betabinom_main","betabinom_2way","betabinom_full3"),
  AIC    = c(AIC(m_bin_full3), AIC(m_bb_main), AIC(m_bb_2way), AIC(m_bb_full3)),
  logLik = c(logLik(m_bin_full3), logLik(m_bb_main), logLik(m_bb_2way), logLik(m_bb_full3)) %>% as.numeric()
)
print(fit_tab)

# Type II tests for the chosen inference model (full 3-way recommended if significant)
cat("\nType II (Beta-binomial full model):\n")
bb_typeII <- car::Anova(m_bb_full3, type = 2)   # Wald χ² per term
print(bb_typeII)

# Residual diagnostics (DHARMa)
set.seed(42)
sim_bb <- simulateResiduals(m_bb_full3, n = 1000)
cat("\nDHARMa tests (beta-binomial full):\n")
print(testUniformity(sim_bb))
print(testDispersion(sim_bb))
# plot(sim_bb)   # uncomment to view plots

# Post-hoc EMMs (response scale) + Tukey; keep or trim as needed
emm_germ <- emmeans(m_bb_full3, ~ Soil * Salinity * Moisture, type = "response")
cat("\nEMMeans (proportion germinated) — first 12 rows:\n")
print(summary(emm_germ, infer = TRUE) %>% as.data.frame() %>% head(12))

cat("\nTukey pairwise (first 25):\n")
print(summary(pairs(emm_germ, adjust = "tukey")) %>% as.data.frame() %>% head(25))

cat("\nCompact Letter Display (CLD):\n")
print(cld(emm_germ, Letters = letters, adjust = "tukey", type = "response"))

#### 3) MEAN GERMINATION TIME (MGT) MODEL --------------------------------------
cat("\n================ MEAN GERMINATION TIME (MGT) MODEL ================\n")

# Build a long-format copy for GerminaR
germsum <- final %>% select(Soil, Salinity, Moisture, Seeds, all_of(day_cols))

# GerminaR summary (per experimental unit aggregated by factors)
gsm <- ger_summary(SeedN = "Seeds", evalName = "D", data = germsum)

# Ensure factors are factors (ger_summary may coerce)
gsm <- gsm %>%
  mutate(
    Soil     = factor(Soil,     levels = levels(final$Soil)),
    Salinity = factor(Salinity, levels = levels(final$Salinity)),
    Moisture = factor(Moisture, levels = levels(final$Moisture))
  ) %>%
  drop_na(mgt)

# Quick sanity prints
cat("\nGerminaR summary columns available:\n")
print(names(gsm))
cat("\nMGT range (days):\n")
print(range(gsm$mgt, na.rm = TRUE))

# Gaussian LM for MGT (transform if needed; start simple)
m_mgt <- lm(mgt ~ Soil * Salinity * Moisture, data = gsm)

# Type II tests (F)
cat("\nType II F-tests (MGT):\n")
print(car::Anova(m_mgt, type = 2))

# Model fit checks for MGT
cat("\nBasic model checks for MGT (normality & heteroscedasticity):\n")
print(performance::check_normality(m_mgt))
print(performance::check_heteroscedasticity(m_mgt))
print(performance::check_collinearity(m_mgt))

# Optional diag plots (uncomment to view)
# par(mfrow = c(2,2)); plot(m_mgt); par(mfrow = c(1,1))

# EMMs + Tukey on MGT
emm_mgt <- emmeans(m_mgt, ~ Soil * Salinity * Moisture)
cat("\nEMMeans (MGT, days) — first 12 rows:\n")
print(summary(emm_mgt, infer = TRUE) %>% as.data.frame() %>% head(12))

cat("\nTukey pairwise for MGT (first 25):\n")
print(summary(pairs(emm_mgt, adjust = "tukey")) %>% as.data.frame() %>% head(25))

cat("\nCLD for MGT:\n")
print(cld(emm_mgt, Letters = letters, adjust = "tukey"))



## ======== MANUSCRIPT SUMMARY (print-only, no files) ==========================
cat("\n================ MANUSCRIPT SUMMARY ==========================\n")

## Germination % (beta-binomial)
cat("\nGermination proportion (beta-binomial GLM):\n")
# Model comparisons you computed
cmp_2vsmain <- anova(m_bb_main, m_bb_2way)
cmp_3vs2    <- anova(m_bb_2way, m_bb_full3)

cat(sprintf("  2-way vs main:  LR Chi^2 = %.2f (df=%d), p = %.3g\n",
            as.numeric(cmp_2vsmain$Chisq[2]),
            as.integer(cmp_2vsmain$`Chi Df`[2]),
            as.numeric(cmp_2vsmain$`Pr(>Chisq)`[2])))

cat(sprintf("  3-way vs 2-way: LR Chi^2 = %.2f (df=%d), p = %.3g\n",
            as.numeric(cmp_3vs2$Chisq[2]),
            as.integer(cmp_3vs2$`Chi Df`[2]),
            as.numeric(cmp_3vs2$`Pr(>Chisq)`[2])))

# AICs
cat("  AICs (lower better):\n")
print(tibble::tibble(
  model  = c("Main","Two-way","Full (3-way)"),
  AIC    = c(AIC(m_bb_main), AIC(m_bb_2way), AIC(m_bb_full3))
))

# Type II table (already computed as bb_typeII)
cat("\n  Type II (Wald Chi^2) – full model terms:\n")
print(bb_typeII)

# A couple of EMMs (extremes)
emm_tab <- summary(emm_germ, infer = TRUE) %>% as.data.frame()
top3 <- emm_tab %>% arrange(desc(prob)) %>% head(3)
bot3 <- emm_tab %>% arrange(prob) %>% head(3)

cat("\n  Highest estimated germination probabilities (top 3):\n")
print(top3[, c("Soil","Salinity","Moisture","prob","asymp.LCL","asymp.UCL")])

cat("\n  Lowest estimated germination probabilities (bottom 3):\n")
print(bot3[, c("Soil","Salinity","Moisture","prob","asymp.LCL","asymp.UCL")])

## MGT (Gaussian LM)
cat("\nMean germination time (Gaussian LM):\n")
a_mgt <- car::Anova(m_mgt, type = 2)
print(a_mgt)

# Fit checks you computed
cat("\n  Fit checks:\n")
print(performance::check_normality(m_mgt))
print(performance::check_heteroscedasticity(m_mgt))

# A few EMMs for context
emm_tab_mgt <- summary(emm_mgt, infer = TRUE) %>% as.data.frame()
cat("\n  EMMs (first 6 rows):\n")
print(head(emm_tab_mgt[, c("Soil","Salinity","Moisture","emmean","lower.CL","upper.CL")]))

cat("\n=============================================================\n")


m_mgt_gamma <- glmmTMB(mgt ~ Soil * Salinity * Moisture,
                       data = gsm, family = Gamma(link = "log"))

set.seed(123)
sim_mgt_g <- simulateResiduals(m_mgt_gamma, n = 1000)
cat("\n[MGT] Gamma GLM DHARMa tests:\n")
print(testUniformity(sim_mgt_g))
print(testDispersion(sim_mgt_g))
# plot(sim_mgt_g)  # optional

cat("\n[MGT] Type II (Gamma GLM):\n")
print(car::Anova(m_mgt_gamma, type = 2))

# Compare Gaussian vs log-Gaussian vs Gamma (pseudo):
cat("\n[MGT] Model comparison (information criteria):\n")
cmp <- tibble::tibble(
  model  = c("Gaussian", "log-Gaussian", "Gamma(log-link)"),
  AIC    = c(AIC(m_mgt), AIC(m_mgt_log), AIC(m_mgt_gamma)),
  logLik = c(logLik(m_mgt), logLik(m_mgt_log), logLik(m_mgt_gamma)) %>% as.numeric()
)
print(cmp)

# EMMs on response scale for Gamma (already on days via log-link)
emm_mgt_gamma <- emmeans(m_mgt_gamma, ~ Soil * Salinity * Moisture, type = "response")
cat("\n[MGT] EMMs (Gamma GLM) — first 6:\n")
print(summary(emm_mgt_gamma, infer = TRUE) %>% as.data.frame() %>% head(6))

#### MGT: add log-Gaussian model + compare all three ###########################

# 1) Fit log-Gaussian (log-transform mgt; you have no zeros, so log() is fine)
gsm <- gsm %>% mutate(mgt_log = log(mgt))

m_mgt_log <- lm(mgt_log ~ Soil * Salinity * Moisture, data = gsm)

cat("\n[MGT] Log-Gaussian (on log-days) — Type II tests:\n")
print(car::Anova(m_mgt_log, type = 2))

cat("\n[MGT] Fit checks (log-Gaussian):\n")
print(performance::check_normality(m_mgt_log))
print(performance::check_heteroscedasticity(m_
                                            
                                            