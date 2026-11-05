
# ============================================================
# SEED SIZE, GERMINATION, AND INUNDATION ANALYSIS
# ============================================================


# ============================================================
# 1. LOAD PACKAGES AND DATA
# ============================================================

library(dplyr)
library(ggplot2)
library(scales)
library(broom)
library(grid)

df <- read.csv("GermSized.csv")


# ============================================================
# 2. CLEAN AND HARMONIZE DATA
# ============================================================

# Harmonize seed source labels
df <- df %>%
  mutate(
    Seed.Source = case_when(
      Seed.Source %in% c("Hester1", "Hester2") ~ "Hester",
      Seed.Source == "Azevedo" ~ "APC",
      TRUE ~ Seed.Source
    )
  )


# ============================================================
# 3. PREPARE VARIABLES
# ============================================================

# Function to convert variables to numeric
to_num <- function(x) {
  suppressWarnings(
    as.numeric(gsub(",", "", trimws(as.character(x))))
  )
}

# Create a unified Size variable
if (!("Size" %in% names(df))) df$Size <- NA_real_

if ("Size" %in% names(df)) {
  df$Size <- ifelse(is.na(df$Size), to_num(df$Size), df$Size)
}

if ("AvgSize" %in% names(df)) {
  df$Size <- ifelse(is.na(df$Size), to_num(df$AvgSize), df$Size)
}

if ("Sized" %in% names(df)) {
  df$Size <- ifelse(is.na(df$Size), to_num(df$Sized), df$Size)
}

# Convert Germ to numeric
if (!("Germ" %in% names(df))) {
  stop("Column 'Germ' not found.")
}

df$Germ <- to_num(df$Germ)

# Convert percentages to proportions if necessary
if (max(df$Germ, na.rm = TRUE) > 1) {
  df$Germ <- df$Germ / 100
}

# Define trial weights
# Prefer Total, then Seeds; otherwise assume 50
if ("Total" %in% names(df)) {
  df$Trials <- to_num(df$Total)
  
} else if ("Seeds" %in% names(df)) {
  df$Trials <- to_num(df$Seeds)
  
} else {
  df$Trials <- 50
}

# Remove incomplete or invalid observations
df <- df %>%
  filter(
    is.finite(Size),
    is.finite(Germ),
    is.finite(Trials),
    Trials > 0
  )


# ============================================================
# 4. CREATE INUNDATION CATEGORIES
# ============================================================

# Frequent inundation:
# SouthMarsh, Kirby, Yampah
#
# Infrequent inundation:
# All remaining sources, including Hester and APC

df <- df %>%
  mutate(
    Inundation = case_when(
      Seed.Source %in% c("SouthMarsh", "Kirby", "Yampah") ~ "Frequent",
      TRUE ~ "Infrequent"
    ),
    Inundation = factor(
      Inundation,
      levels = c("Infrequent", "Frequent")
    )
  )

# Set plot colors
COL_INFREQ <- "#1b9e77"
COL_FREQ   <- "#d95f02"

# Establish common x-axis range
x_min <- min(df$Size, na.rm = TRUE)
x_max <- max(df$Size, na.rm = TRUE)
x_span <- x_max - x_min


# ============================================================
# 5. EXPLORATORY DATA ANALYSIS
# ============================================================

# Number of observations per seed source
df %>%
  count(Seed.Source, Inundation)

# Number of unique seed sources per inundation category
df %>%
  distinct(Seed.Source, Inundation) %>%
  count(Inundation)

# Summary statistics by seed source
source_summary <- df %>%
  group_by(Seed.Source, Inundation) %>%
  summarise(
    Mean_Germ = mean(Germ),
    SD_Germ = sd(Germ),
    N = n(),
    .groups = "drop"
  )

print(source_summary)

# Overall germination summary by inundation category
inundation_summary <- df %>%
  group_by(Inundation) %>%
  summarise(
    Mean_Germ = mean(Germ),
    SD_Germ = sd(Germ),
    N = n(),
    .groups = "drop"
  )

print(inundation_summary)


# ============================================================
# 6. OVERALL WEIGHTED LINEAR REGRESSION
# ============================================================

# Fit weighted linear regression
fit_all <- lm(
  Germ ~ Size,
  data = df,
  weights = Trials
)

# Model summary
summary(fit_all)

# Extract regression coefficients and weighted R-squared
s_all  <- summary(fit_all)
b0_all <- unname(coef(fit_all)[1])
b1_all <- unname(coef(fit_all)[2])
r2_all <- s_all$r.squared

# Create equation label
lab_all <- data.frame(
  x = x_max - 0.05 * x_span,
  y = 0.95,
  label = paste0(
    "y = ", sprintf("%.3f", b0_all),
    ifelse(b1_all >= 0, " + ", " - "),
    sprintf("%.3f", abs(b1_all)), "·x",
    "\nR² (weighted) = ", sprintf("%.2f", r2_all)
  )
)


# ============================================================
# 7. INUNDATION-SPECIFIC WEIGHTED LINEAR REGRESSIONS
# ============================================================

# Fit separate weighted regressions for each inundation category
fit_by <- df %>%
  filter(Inundation %in% c("Infrequent", "Frequent")) %>%
  group_by(Inundation) %>%
  do({
    fit <- lm(
      Germ ~ Size,
      data = .,
      weights = Trials
    )
    
    s <- summary(fit)
    
    tibble(
      b0 = unname(coef(fit)[1]),
      b1 = unname(coef(fit)[2]),
      r2 = s$r.squared
    )
  }) %>%
  ungroup() %>%
  mutate(
    label = paste0(
      "y = ", sprintf("%.3f", b0),
      ifelse(b1 >= 0, " + ", " - "),
      sprintf("%.3f", abs(b1)), "·x",
      "\nR² (weighted) = ", sprintf("%.2f", r2)
    )
  )

# Set equation label positions
label_pos <- data.frame(
  Inundation = c("Infrequent", "Frequent"),
  x = x_max - 0.05 * x_span,
  y = c(0.95, 0.85)
)

fit_by <- left_join(
  fit_by,
  label_pos,
  by = "Inundation"
)

# View regression results
print(fit_by)


# ============================================================
# 8. PLOT A: OVERALL WEIGHTED LINEAR REGRESSION
# ============================================================

LeeSupplementalFigure5 <- ggplot(
  df,
  aes(x = Size, y = Germ, color = Inundation)
) +
  geom_point(
    size = 2.6,
    alpha = 0.95
  ) +
  geom_smooth(
    aes(weight = Trials),
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "grey80",
    linewidth = 0.9
  ) +
  geom_label(
    data = lab_all,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = label),
    fill = "white",
    label.size = 0.25,
    size = 3.6,
    hjust = 1
  ) +
  scale_color_manual(
    values = c(
      "Infrequent" = COL_INFREQ,
      "Frequent" = COL_FREQ
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = "Seed Size (mm)",
    y = "Germination Percentage (%)",
    color = "Inundation"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

# Display plot
print(LeeSupplementalFigure5)

# Save plot
ggsave(
  "LeeSupplementalFigure5.tiff",
  plot = LeeSupplementalFigure5,
  width = 5,
  height = 5,
  dpi = 300,
  units = "in"
)


# ============================================================
# 9. PLOT B: INUNDATION-SPECIFIC WEIGHTED REGRESSIONS
# ============================================================

p_bybin <- ggplot(
  df,
  aes(x = Size, y = Germ, color = Inundation)
) +
  geom_point(
    size = 2.6,
    alpha = 0.95
  ) +
  geom_smooth(
    aes(weight = Trials),
    method = "lm",
    se = TRUE,
    linewidth = 0.9
  ) +
  geom_label(
    data = fit_by,
    inherit.aes = FALSE,
    aes(
      x = x,
      y = y,
      label = label,
      color = Inundation
    ),
    fill = "white",
    label.size = 0.25,
    size = 3.6,
    show.legend = FALSE,
    hjust = 1
  ) +
  scale_color_manual(
    values = c(
      "Infrequent" = COL_INFREQ,
      "Frequent" = COL_FREQ
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    title = "Seed Size vs Germination by Inundation (Weighted Linear Fits)",
    x = "Seed Size (mm)",
    y = "Germination (proportion)",
    color = "Inundation"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

# Display plot
print(p_bybin)


# ============================================================
# 10. WEIGHTED LINEAR MODEL DIAGNOSTICS
# ============================================================

# Diagnostic plots for overall regression
par(mfrow = c(2, 2))
plot(fit_all)
par(mfrow = c(1, 1))

# Overall model summary
summary(fit_all)


# ============================================================
# 11. QUASIBINOMIAL MODEL: INUNDATION EFFECT
# ============================================================

# Test whether germination differs between inundation categories
mod_inundation <- glm(
  Germ ~ Inundation,
  family = quasibinomial,
  data = df
)

# Model summary
summary(mod_inundation)

# Check overdispersion
deviance(mod_inundation) / df.residual(mod_inundation)

# Calculate odds ratios and confidence intervals
exp(
  cbind(
    OR = coef(mod_inundation),
    confint.default(mod_inundation)
  )
)


# ============================================================
# 12. SOURCE-LEVEL GERMINATION COMPARISON
# ============================================================

# Calculate average germination for each source
source_means <- df %>%
  group_by(Seed.Source, Inundation) %>%
  summarise(
    Mean_Germ = mean(Germ),
    .groups = "drop"
  )

# Compare source-level mean germination between inundation groups
t.test(
  Mean_Germ ~ Inundation,
  data = source_means
)
