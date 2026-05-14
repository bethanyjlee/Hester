# ===============================
# Packages
# ===============================
library(dplyr)
library(ggplot2)
library(scales)
library(broom)
library(grid)   # for unit()

# ===============================
# 1) Load & harmonize labels
# ===============================
df <- read.csv("GermSized.csv")

# Merge Hester1/Hester2 -> Hester; rename Azevedo -> APC
df <- df %>%
  mutate(
    Seed.Source = case_when(
      Seed.Source %in% c("Hester1","Hester2") ~ "Hester",
      Seed.Source == "Azevedo"                ~ "APC",
      TRUE                                    ~ Seed.Source
    )
  )

# ===============================
# 2) Build robust Size, Germ (0–1), Trials (weights)
# ===============================
to_num <- function(x) suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(x)))))

# Create a single 'Size' from any of: Size, AvgSize, Sized
if (!("Size" %in% names(df))) df$Size <- NA_real_
if ("Size"    %in% names(df)) df$Size <- ifelse(is.na(df$Size), to_num(df$Size), df$Size)
if ("AvgSize" %in% names(df)) df$Size <- ifelse(is.na(df$Size), to_num(df$AvgSize), df$Size)
if ("Sized"   %in% names(df)) df$Size <- ifelse(is.na(df$Size), to_num(df$Sized),   df$Size)

# Germ to 0–1
if (!("Germ" %in% names(df))) stop("Column 'Germ' not found.")
df$Germ <- to_num(df$Germ)
if (is.finite(suppressWarnings(max(df$Germ, na.rm = TRUE))) && max(df$Germ, na.rm = TRUE) > 1) {
  df$Germ <- df$Germ / 100
}

# Trials (weights): prefer 'Total' or 'Seeds'; else assume 50
if ("Total" %in% names(df)) {
  df$Trials <- to_num(df$Total)
} else if ("Seeds" %in% names(df)) {
  df$Trials <- to_num(df$Seeds)
} else {
  df$Trials <- 50
}

# Keep only complete rows
df <- df %>%
  filter(is.finite(Size), is.finite(Germ), is.finite(Trials), Trials > 0)

# ===============================
# 3) Inundation bins only (no site factors in plots)
# Frequent: SouthMarsh, Kirby, Yampah; others (incl. Hester/APC) Infrequent
# ===============================
df <- df %>%
  mutate(
    Inundation = case_when(
      Seed.Source %in% c("SouthMarsh","Kirby","Yampah") ~ "Frequent",
      TRUE                                              ~ "Infrequent"
    )
  )

# Colors to match earlier figs
COL_INFREQ <- "#1b9e77"
COL_FREQ   <- "#d95f02"

# Common x-range helpers
x_min <- min(df$Size, na.rm = TRUE); x_max <- max(df$Size, na.rm = TRUE)
x_span <- x_max - x_min

# ===============================
# 4) PLOT A — Overall weighted linear regression
# One line with CI; points colored by Inundation; show equation & weighted R²
# ===============================
fit_all <- lm(Germ ~ Size, data = df, weights = Trials)
s_all   <- summary(fit_all)
b0_all  <- unname(coef(fit_all)[1])
b1_all  <- unname(coef(fit_all)[2])
r2_all  <- s_all$r.squared   # weighted R²

lab_all <- data.frame(
  x = x_max - 0.05 * x_span,   # near the right edge
  y = 0.95,
  label = paste0(
    "y = ", sprintf("%.3f", b0_all),
    ifelse(b1_all >= 0, " + ", " - "), sprintf("%.3f", abs(b1_all)), "·x",
    "\nR² (weighted) = ", sprintf("%.2f", r2_all)
  )
)

LeeSupplementalFigure5 <- ggplot(df, aes(x = Size, y = Germ, color = Inundation)) +
  geom_point(size = 2.6, alpha = 0.95) +
  geom_smooth(aes(weight = Trials), method = "lm", se = TRUE, color = "black", fill = "grey80", size = 0.9) +
  geom_label(
    data = lab_all, inherit.aes = FALSE,
    aes(x = x, y = y, label = label),
    fill = "white", label.size = 0.25, size = 3.6,
    hjust = 1        # right-align text box
  ) +
  scale_color_manual(values = c("Infrequent" = COL_INFREQ, "Frequent" = COL_FREQ)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = "Seed Size (mm)", y = "Germination Percentage (%)", color = "Inundation") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

print(p_overall)

ggsave("LeeSupplementalFigure5.tiff", width =5, height = 5, dpi = 300, units = "in")

# ===============================
# 5) PLOT B — Per-inundation weighted linear regressions
# Two lines (one per bin) with their own CIs + equations & weighted R²s
# ===============================
# Per-bin weighted linear models
fit_by <- df %>%
  filter(Inundation %in% c("Infrequent","Frequent")) %>%
  group_by(Inundation) %>%
  do({
    fit <- lm(Germ ~ Size, data = ., weights = Trials)
    s   <- summary(fit)
    tibble(
      b0 = unname(coef(fit)[1]),
      b1 = unname(coef(fit)[2]),
      r2 = s$r.squared
    )
  }) %>%
  ungroup() %>%
  mutate(
    label = paste0("y = ", sprintf("%.3f", b0),
                   ifelse(b1 >= 0, " + ", " - "), sprintf("%.3f", abs(b1)), "·x",
                   "\nR² (weighted) = ", sprintf("%.2f", r2))
  )

label_pos <- data.frame(
  Inundation = c("Infrequent","Frequent"),
  x = x_max - 0.05 * x_span,   # near the right edge
  y = c(0.95, 0.85)
)
fit_by <- left_join(fit_by, label_pos, by = "Inundation")

p_bybin <- ggplot(df, aes(x = Size, y = Germ, color = Inundation)) +
  geom_point(size = 2.6, alpha = 0.95) +
  geom_smooth(aes(weight = Trials), method = "lm", se = TRUE, size = 0.9) +
  geom_label(
    data = fit_by, inherit.aes = FALSE,
    aes(x = x, y = y, label = label, color = Inundation),
    fill = "white", label.size = 0.25, size = 3.6, show.legend = FALSE,
    hjust = 1        # right-align text boxes
  ) +
  scale_color_manual(values = c("Infrequent" = COL_INFREQ, "Frequent" = COL_FREQ)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(title = "Seed Size vs Germination by Inundation (Weighted Linear Fits)",
       x = "Seed Size (mm)", y = "Germination (proportion)", color = "Inundation") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

print(p_bybin)
