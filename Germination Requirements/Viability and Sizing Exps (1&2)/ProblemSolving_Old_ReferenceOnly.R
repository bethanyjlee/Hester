# ---- Libraries ----
library(tidyverse)
library(pscl)     # pR2() for pseudo-R²
library(broom)    # tidy outputs (optional)
# library(car)    # if you want Levene's test later

# ---- Load data ----
sized <- read.csv("GermSized.csv")

# ---- Auto-detect key columns (edit if needed) ----
cn <- names(sized)

site_var <- cn[grepl("seed.source|^site$|source", cn, ignore.case = TRUE)][1]
size_var <- cn[grepl("^size$|avgsize|sized|seed.size", cn, ignore.case = TRUE)][1]
germ_var <- cn[grepl("^germ$|germination|perc", cn, ignore.case = TRUE)][1]
total_var <- cn[grepl("^total$|total.seed|n.seed|seeds.total", cn, ignore.case = TRUE)][1]

stopifnot(!is.na(site_var), !is.na(size_var), !is.na(germ_var))  # fail early if not found

dat <- sized %>%
  rename(
    Site = all_of(site_var),
    Size = all_of(size_var),
    Germ_raw = all_of(germ_var)
  ) %>%
  mutate(
    # If Germ is in 0–100, convert to 0–1
    Germ = ifelse(Germ_raw > 1, Germ_raw / 100, Germ_raw),
    Total = if (!is.na(total_var)) .data[[total_var]] else 50
  ) %>%
  filter(!is.na(Site), !is.na(Size), !is.na(Germ), !is.na(Total))

# ---- 1) Size differences among sites ----
aov_size <- aov(Size ~ Site, data = dat)
cat("\n=== Size ~ Site ANOVA ===\n")
print(summary(aov_size))
cat("\n--- Tukey post-hoc ---\n")
print(TukeyHSD(aov_size))

# (Optional nonparametric fallback if needed)
kruskal.test(Size ~ Site, data = dat)

# ---- 2) GLMs for germination ~ size (with/without site) ----
# Build successes/failures from proportion + totals
dat <- dat %>%
  mutate(
    Success = round(Germ * Total),
    Failure = round((1 - Germ) * Total)
  )

glm_site   <- glm(cbind(Success, Failure) ~ Size + Site, data = dat, family = binomial)
glm_nosite <- glm(cbind(Success, Failure) ~ Size,        data = dat, family = binomial)

cat("\n=== GLM with Site ===\n");   print(summary(glm_site))
cat("\n=== GLM without Site ===\n"); print(summary(glm_nosite))

cat("\n=== Model comparison (AIC) ===\n")
print(AIC(glm_site, glm_nosite))

cat("\n=== Likelihood-ratio test (does Site improve fit?) ===\n")
print(anova(glm_nosite, glm_site, test = "Chisq"))

## --- Pseudo-R²s computed manually (robust with cbind(...) syntax) ---

# Fit explicit null models (intercept only) using the same response
glm_null       <- glm(cbind(Success, Failure) ~ 1, data = dat, family = binomial)
glm_null_nosite <- glm_null  # same null for both comparisons

# Log-likelihoods
ll_null       <- as.numeric(logLik(glm_null))
ll_site       <- as.numeric(logLik(glm_site))
ll_nosite     <- as.numeric(logLik(glm_nosite))

# Number of estimated parameters (including intercept)
k_site        <- length(coef(glm_site))
k_nosite      <- length(coef(glm_nosite))

# McFadden R²
r2_mcf_site   <- 1 - (ll_site   / ll_null)
r2_mcf_nosite <- 1 - (ll_nosite / ll_null)

# Adjusted McFadden R² (penalizes extra parameters)
r2_mcf_adj_site   <- 1 - ((ll_site   - k_site)   / ll_null)
r2_mcf_adj_nosite <- 1 - ((ll_nosite - k_nosite) / ll_null)

# Print them
cbind(
  model = c("with site", "no site"),
  McFadden = round(c(r2_mcf_site, r2_mcf_nosite), 3),
  McFadden_adj = round(c(r2_mcf_adj_site, r2_mcf_adj_nosite), 3)
)
# Nice positions for labels
x_left  <- quantile(dat$Size, 0.02, na.rm = TRUE)
y_top   <- 0.98
gap     <- 0.08

ggplot(dat, aes(x = Size, y = Germ, color = Site)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Seed Size (mm)",
    y = "Germination (%)",
    color = "Site"
  ) +
  annotate("text", x = x_left, y = y_top,
           label = paste0("McFadden R² = ", round(r2_mcf_site, 3)),
           hjust = 0, vjust = 1) +
  annotate("text", x = x_left, y = y_top - gap,
           label = paste0("Adj. McFadden R² = ", round(r2_mcf_adj_site, 3)),
           hjust = 0, vjust = 1) +
  theme_classic()

#### checking binomial distribution 

family(glm_site)
# Residual deviance / residual degrees of freedom
dispersion <- sum(residuals(glm_site, type = "deviance")^2) / df.residual(glm_site)
dispersion
par(mfrow = c(2,2))
plot(glm_site)  # standard diagnostic plots
library(DHARMa)
sim_res <- simulateResiduals(glm_site, n = 1000)
plot(sim_res)
glm_site_quasi <- glm(cbind(Success, Failure) ~ Size + Site, data = dat, family = quasibinomial)
summary(glm_site_quasi)
AIC(glm_site, glm_site_quasi)

