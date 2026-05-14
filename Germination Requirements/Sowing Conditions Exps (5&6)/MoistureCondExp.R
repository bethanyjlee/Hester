#### Packages ----
library(dplyr)
library(car)        # Anova()
library(emmeans)    # emmeans, pairs
library(multcomp)   # cld() for compact letters
library(broom)
library(ggplot2)

setwd("~/R data/Chp1PWSeeds/Zea and Wesley Exp")

#### Data & prep ----


dat <- read.csv("Zeanew.csv")

# If any 'air' rows remain, drop them (safe even if none are present)
dat <- dat %>% filter(is.na(exposure) | exposure == "container")

dat <- dat %>%
  mutate(
    # Use the levels present in the file
    Moisture_level_QL = factor(Moisture_level_QL,
                               levels = c("Lo","MdLo","MdHi","Hi","Water")),
    Photo_period      = factor(Photo_period)
  ) %>%
  # keep rows with complete info needed for the model
  filter(!is.na(MaxGerm), !is.na(Total),
         !is.na(Moisture_level_QL), !is.na(Photo_period), Total > 0) %>%
  droplevels()

# Quick sanity check: both factors must have >= 2 levels
stopifnot(nlevels(dat$Moisture_level_QL) >= 2,
          nlevels(dat$Photo_period)      >= 2)

#### Main-effects model ----
m_main <- glm(cbind(MaxGerm, Total - MaxGerm) ~ Moisture_level_QL + Photo_period,
              family = binomial(link = "logit"), data = dat)

# Overdispersion check
disp <- sum(residuals(m_main, type="deviance")^2) / df.residual(m_main)

if (disp > 1.5) {
  message(sprintf("Overdispersion detected (%.2f) -> using quasibinomial.", disp))
  m_final <- glm(cbind(MaxGerm, Total - MaxGerm) ~ Moisture_level_QL + Photo_period,
                 family = quasibinomial(link = "logit"), data = dat)
  fam_used <- "quasibinomial"
} else {
  # If dispersion is tame, optionally test interaction via AIC
  m_int <- update(m_main, . ~ Moisture_level_QL * Photo_period)
  aics <- AIC(m_main, m_int)
  keep_int <- (aics$AIC[2] + 2) < aics$AIC[1]  # keep if improves by >2 AIC units
  m_final <- if (keep_int) m_int else m_main
  fam_used <- "binomial"
  message(if (keep_int) "Keeping interaction (AIC improved)." else "Using main-effects (no AIC improvement).")
}

#### Global tests ----
tests <- car::Anova(m_final, type = 2)  # Type II; for quasi, this is Wald w/ dispersion
print(tests)

# Identify significant terms
sig_terms <- rownames(tests)[which(tests$`Pr(>Chisq)` < 0.05)]
cat("\nSignificant terms (α=0.05): ",
    ifelse(length(sig_terms)==0, "None", paste(sig_terms, collapse=", ")), "\n")

#### EMMs & Tukey only for significant factors ----
emm_out <- list(); pairs_out <- list(); cld_out <- list()

if ("Moisture_level_QL" %in% sig_terms) {
  emm_out$moist <- emmeans(m_final, ~ Moisture_level_QL, type = "response")
  pairs_out$moist <- pairs(emm_out$moist, adjust = "tukey")
  cld_out$moist <- multcomp::cld(emm_out$moist, Letters = letters,
                                 adjust = "tukey", type = "response")
  cat("\nMoisture emmeans (prob, 95% CI) + Tukey groups:\n")
  print(as.data.frame(cld_out$moist)[, c("Moisture_level_QL","prob","asymp.LCL","asymp.UCL",".group")])
  cat("\nMoisture Tukey pairs:\n"); print(pairs_out$moist)
}

if ("Photo_period" %in% sig_terms) {
  emm_out$photo <- emmeans(m_final, ~ Photo_period, type = "response")
  pairs_out$photo <- pairs(emm_out$photo, adjust = "tukey")
  cat("\nPhotoperiod emmeans (prob, 95% CI):\n")
  print(as.data.frame(emm_out$photo)[, c("Photo_period","prob","asymp.LCL","asymp.UCL")])
  cat("\nPhotoperiod Tukey pair:\n"); print(pairs_out$photo)
}

#### (Optional) Nice compact tables for Results text ----
overall_tbl <- tests %>%
  tibble::rownames_to_column("Term") %>%
  transmute(Term,
            Stat = round(Chisq, 2),
            df = Df,
            p = signif(`Pr(>Chisq)`, 3))
print(overall_tbl)



#### Test for removing Photoperiod #####
# ---- Packages ----
library(dplyr)
library(car)        # Anova()
library(emmeans)    # (not strictly needed here, but useful if you want posthocs)
# Optional for overdispersion-aware simulation; we handle fallback if not present:
has_extrad <- requireNamespace("extraDistr", quietly = TRUE)

# ---- Data prep (matches your successful run) ----
dat <- read.csv("Zeanew.csv") %>%
  filter(is.na(exposure) | exposure == "container") %>%     # ensure no 'air'
  mutate(
    Moisture_level_QL = factor(Moisture_level_QL,
                               levels = c("Lo","MdLo","MdHi","Hi","Water")),
    Photo_period      = factor(Photo_period)
  ) %>%
  filter(!is.na(MaxGerm), !is.na(Total),
         !is.na(Moisture_level_QL), !is.na(Photo_period), Total > 0) %>%
  droplevels()

stopifnot(nlevels(dat$Moisture_level_QL) >= 2,
          nlevels(dat$Photo_period)      >= 2)

# ---- Fit the two candidate models on the real data ----
m_phot <- glm(cbind(MaxGerm, Total - MaxGerm) ~ Moisture_level_QL + Photo_period,
              family = binomial(link = "logit"), data = dat)
# Detect overdispersion (you observed ~1.74 earlier)
disp <- sum(residuals(m_phot, type="deviance")^2) / df.residual(m_phot)

# We use the *no-photoperiod* fit as the data-generating truth (photoperiod ~ null)
m_nophot <- glm(cbind(MaxGerm, Total - MaxGerm) ~ Moisture_level_QL,
                family = binomial(link = "logit"), data = dat)

# Linear predictor & probabilities under the null (no photoperiod effect)
lp_null <- model.matrix(~ Moisture_level_QL, dat) %*% coef(m_nophot)
p_null  <- plogis(as.numeric(lp_null))

# ---- Convert dispersion to a beta-binomial rho (if we can) ----
# Var_BB = n*p*(1-p) * (1 + (n-1)*rho)
# Quasi says Var_quasi = phi * n*p*(1-p). Equate at typical n_bar: rho ≈ (phi - 1)/(n_bar - 1)
n_bar <- mean(dat$Total)
phi   <- max(disp, 1)   # at least 1
rho   <- max(min((phi - 1) / max(n_bar - 1, 1e-9), 0.99), 0)  # clamp to [0, 0.99]

message(sprintf("Dispersion = %.2f; using rho ≈ %.3f for beta-binomial (if available).", disp, rho))

# ---- Simulation function ----
sim_power <- function(B = 2000, alpha = 0.05, use_beta_binomial = TRUE) {
  # storage
  reject_moist_withP  <- logical(B)
  reject_moist_noP    <- logical(B)
  reject_photo        <- logical(B)  # Type I for photoperiod when truth is null
  
  for (b in seq_len(B)) {
    # Generate counts under the null photoperiod effect, moisture means from m_nophot
    if (use_beta_binomial && has_extrad && rho > 0) {
      # extraDistr::rbbinom uses parameter rho
      y <- extraDistr::rbbinom(n = nrow(dat), size = dat$Total, prob = p_null, rho = rho)
    } else {
      # binomial fallback (conservative variance if phi>1)
      y <- rbinom(n = nrow(dat), size = dat$Total, prob = p_null)
    }
    
    sim <- dat
    sim$MaxGerm <- y
    
    # Fit + test: with Photoperiod
    fit_with <- glm(cbind(MaxGerm, Total - MaxGerm) ~ Moisture_level_QL + Photo_period,
                    family = binomial(link = "logit"), data = sim)
    # Use quasi if original showed notable overdispersion
    if (disp > 1.5) {
      fit_with <- glm(cbind(MaxGerm, Total - MaxGerm) ~ Moisture_level_QL + Photo_period,
                      family = quasibinomial(link = "logit"), data = sim)
    }
    a_with <- car::Anova(fit_with, type = 2)
    p_moist_with <- a_with[ "Moisture_level_QL", grep("^Pr\\(>", colnames(a_with)) ]
    p_photo      <- a_with[ "Photo_period",       grep("^Pr\\(>", colnames(a_with)) ]
    reject_moist_withP[b] <- (p_moist_with < alpha)
    reject_photo[b]       <- (p_photo       < alpha)
    
    # Fit + test: without Photoperiod
    fit_no <- glm(cbind(MaxGerm, Total - MaxGerm) ~ Moisture_level_QL,
                  family = binomial(link = "logit"), data = sim)
    if (disp > 1.5) {
      fit_no <- glm(cbind(MaxGerm, Total - MaxGerm) ~ Moisture_level_QL,
                    family = quasibinomial(link = "logit"), data = sim)
    }
    a_no <- car::Anova(fit_no, type = 2)
    p_moist_no <- a_no[ "Moisture_level_QL", grep("^Pr\\(>", colnames(a_no)) ]
    reject_moist_noP[b] <- (p_moist_no < alpha)
  }
  
  data.frame(
    B = B,
    alpha = alpha,
    dispersion_used = disp,
    rho_used        = if (has_extrad) rho else 0,
    generator       = if (has_extrad && rho > 0) "beta-binomial" else "binomial",
    power_moist_with_photoperiod = mean(reject_moist_withP, na.rm = TRUE),
    power_moist_no_photoperiod   = mean(reject_moist_noP,   na.rm = TRUE),
    type1_photoperiod_when_null  = mean(reject_photo,       na.rm = TRUE)
  )
}

# ---- Run it ----
set.seed(123)
res <- sim_power(B = 2000, alpha = 0.05, use_beta_binomial = TRUE)
print(res)

# Quick delta in power for moisture:
delta <- with(res, power_moist_no_photoperiod - power_moist_with_photoperiod)
message(sprintf("Δ Power (Moisture; no-Photoperiod minus with-Photoperiod) = %.3f", delta))


##### statistical summary table 
#### ================================
#### Moisture-only model + per-level table
#### ================================

# Refit without Photoperiod
m_moist <- glm(cbind(MaxGerm, Total - MaxGerm) ~ Moisture_level_QL,
               family = binomial(link = "logit"), data = dat)

# Overdispersion check (use same rule of thumb as above)
disp_m <- sum(residuals(m_moist, type="deviance")^2) / df.residual(m_moist)
if (disp_m > 1.5) {
  message(sprintf("Overdispersion for moisture-only (%.2f) -> using quasibinomial.", disp_m))
  m_moist <- glm(cbind(MaxGerm, Total - MaxGerm) ~ Moisture_level_QL,
                 family = quasibinomial(link = "logit"), data = dat)
  fam_used_moist <- "quasibinomial"
} else {
  fam_used_moist <- "binomial"
}

# Global (Type II Wald) test for Moisture
tests_moist <- car::Anova(m_moist, type = 2)
print(tests_moist)

## --- Pull χ² and p for caption/text (note: "LR Chisq" name) ---
chi2_m   <- tests_moist[ "Moisture_level_QL", "LR Chisq" ]
p_main_m <- tests_moist[ "Moisture_level_QL", "Pr(>Chisq)" ]

## --- EMMs on link (logit) and response (prob) scales ---
# Use summary(..., infer=TRUE) so we get CIs; GLMs give asymptotic CIs => asymp.LCL/UCL
emm_link_m <- emmeans::emmeans(m_moist, ~ Moisture_level_QL, type = "link") |>
  summary(infer = TRUE) |>
  as.data.frame() |>
  dplyr::rename(
    Level     = Moisture_level_QL,
    logit_est = emmean,
    logit_SE  = SE,
    logit_LCL = asymp.LCL,
    logit_UCL = asymp.UCL
  )

emm_resp_m <- emmeans::emmeans(m_moist, ~ Moisture_level_QL, type = "response") |>
  summary(infer = TRUE) |>
  as.data.frame() |>
  dplyr::rename(
    Level      = Moisture_level_QL,
    model_prob = prob,
    model_SE   = SE,
    model_LCL  = asymp.LCL,
    model_UCL  = asymp.UCL
  ) |>
  dplyr::select(Level, model_prob, model_LCL, model_UCL)  # explicit dplyr::

## --- Optional: compact letters for figure captions (sidak is fine for >2 groups) ---
cld_m <- tryCatch({
  multcomp::cld(emmeans::emmeans(m_moist, ~ Moisture_level_QL, type = "response"),
                Letters = letters, adjust = "tukey", type = "response") |>
    as.data.frame() |>
    dplyr::select(Moisture_level_QL, .group) |>
    dplyr::rename(Level = Moisture_level_QL, Tukey_group = .group)
}, error = function(e) NULL)

## --- Descriptives per level ---
desc_m <- dat |>
  dplyr::group_by(Moisture_level_QL) |>
  dplyr::summarise(
    n_cups    = dplyr::n(),
    seeds_sum = sum(Total,    na.rm = TRUE),
    germ_sum  = sum(MaxGerm,  na.rm = TRUE),
    obs_mean  = ifelse(seeds_sum > 0, germ_sum / seeds_sum, NA_real_),
    .groups = "drop"
  ) |>
  dplyr::rename(Level = Moisture_level_QL)

## --- (Optional) VWC summary if present ---
if ("VWC" %in% names(dat)) {
  vwc_summ <- dat |>
    dplyr::group_by(Moisture_level_QL) |>
    dplyr::summarise(
      VWC_mean = mean(VWC, na.rm = TRUE),
      VWC_sd   = sd(VWC,   na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::rename(Level = Moisture_level_QL)
} else vwc_summ <- NULL

## --- Merge to manuscript-ready table ---
tab_m <- emm_link_m |>
  dplyr::left_join(emm_resp_m, by = "Level") |>
  dplyr::left_join(desc_m,     by = "Level") |>
  { if (!is.null(cld_m)) dplyr::left_join(., cld_m, by = "Level") else . } |>
  { if (!is.null(vwc_summ)) dplyr::left_join(., vwc_summ, by = "Level") else . } |>
  dplyr::mutate(
    `Observed mean`    = sprintf("%.1f%%", 100*obs_mean),
    `EMM (prob)`       = sprintf("%.1f%%", 100*model_prob),
    `95%% CI (prob)`   = paste0(sprintf("%.1f%%", 100*model_LCL), " – ",
                                sprintf("%.1f%%", 100*model_UCL)),
    `Estimate (logit)` = sprintf("%.2f", logit_est),
    `SE (logit)`       = sprintf("%.2f", logit_SE),
    `95%% CI (logit)`  = paste0("[", sprintf("%.2f", logit_LCL), ", ",
                                sprintf("%.2f", logit_UCL), "]")
  ) |>
  dplyr::select(
    Level,
    n_cups, seeds_sum,
    dplyr::any_of(c("VWC_mean","VWC_sd")),
    `Observed mean`, `EMM (prob)`, `95% CI (prob)`,
    `Estimate (logit)`, `SE (logit)`, `95% CI (logit)`,
    dplyr::any_of("Tukey_group")
  ) |>
  dplyr::rename(
    Treatment     = Level,
    `n (cups)`    = n_cups,
    `Total seeds` = seeds_sum,
    `VWC mean`    = VWC_mean,
    `VWC sd`      = VWC_sd
  ) |>
  dplyr::arrange(factor(Treatment, levels = c("Lo","MdLo","MdHi","Hi","Water")))

cat("\n=== Moisture-only per-treatment table ===\n")
print(tab_m, row.names = FALSE)

utils::write.csv(tab_m, "Table_Moisture_Per_Treatment.csv", row.names = FALSE)

cat(sprintf(
  "\nMain effect of Moisture (family=%s): LR chi^2 = %.2f, p = %s; dispersion = %.2f\n",
  fam_used_moist, chi2_m, format.pval(p_main_m, digits = 3), disp_m
))


moist_levels <- c("Lo", "MdLo", "MdHi", "Hi", "Water",
                 setdiff(sort(unique(dat$Soil)), c("Low","Medium-Low","Medium-High","High","Water")))


dat <- dat %>% mutate(Moisture_level = factor(Moisture_level_QL, levels = moist_levels))

#### make supplemental plot #####

waterplot<- ggplot(dat, aes(x=Moisture_level_QL, y=MaxGerm, fill =Treatment)) +
       geom_boxplot(outlier.shape = 21, alpha = 0.9) +
  scale_x_discrete(labels = c(Lo    = "Low", MdLo  = "Medium-Low",
                              MdHi  = "Medium-High", Hi  = "High",
                              Water = "Water")) +
  scale_fill_manual(values = c(
    "Moisture Gradient" = "#6C757D",  
    "Water"             = "#2C7FB8")) +
         labs(x = "Saturation Treatment", y = "Germination (%)") +
         scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
         theme_bw() +
         theme(
           panel.grid.minor = element_blank(),
           plot.title = element_text(face = "bold"),
           legend.position = "none",
           strip.text = element_text(face = "bold"),
           axis.title.x = element_text(size = 12),
           axis.title.y = element_text(size = 12))


ggsave("LeeFigure4.tiff", width = 5, height = 4, dpi = 300, units = "in")
