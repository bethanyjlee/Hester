#### ChatGPT Recommendations
####prompt
####fromGPTOLDREC 9/11

# ============================================================
# Mega2New.csv (DAILY NEW germination) -> GerminaR, GLMs, viridis boxplots
# ============================================================

# ---- Packages ----
need_pkgs <- c("GerminaR","dplyr","ggplot2","readr","stringr")
for (p in need_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  library(GerminaR)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

# ---- Settings ----
input_path <- "Mega2New.csv"     # change if needed
out_dir    <- "plots"; dir.create(out_dir, showWarnings = FALSE)

# Axis/factor orders (edit if your labels differ)
moisture_levels <- c("Low","Medium","High")
salinity_levels <- c("Freshwater","Saltwater")
soil_levels     <- c("HesterMimic","PottingSoil")   # <-- Hester first

soil_labels <- c(HesterMimic = "Hester Mimic",
                 PottingSoil = "Potting Soil")


# ---- 1) Read & standardize columns ----
raw <- read_csv(input_path, show_col_types = FALSE)

# Trim whitespace from ALL column names (fixes 'Salinity ' -> 'Salinity')
names(raw) <- trimws(names(raw))

# Ensure mandatory columns exist
needed <- c("Soil","Salinity","Moisture")
missing_needed <- setdiff(needed, names(raw))
if (length(missing_needed)) {
  stop("Missing required columns in Mega2New.csv: ", paste(missing_needed, collapse=", "))
}

# Seeds: use provided column if present, else default to 146
if (!("Seeds" %in% names(raw))) raw$Seeds <- 146
raw$Seeds <- suppressWarnings(as.numeric(raw$Seeds))
raw$Seeds[is.na(raw$Seeds)] <- 146

# ---- 2) Identify day columns (D5, D7, ... are OK) and coerce numeric ----
day_cols <- grep("^D\\d+$", names(raw), value = TRUE)
if (length(day_cols) == 0) stop("No day columns like D5, D7, ... found in Mega2New.csv")

to_num <- function(x) suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(x)), fixed = TRUE)))

raw <- raw %>%
  mutate(
    across(all_of(day_cols), ~ pmax(to_num(.), 0)),  # daily NEW counts; non-negative
    Seeds = as.numeric(Seeds)
  )

# ---- 3) GerminaR summary (daily NEW -> grs & mgt) ----
# evalName is the prefix of day columns; ours are D5..D23 -> "D"
gsm <- GerminaR::ger_summary(
  SeedN    = "Seeds",
  evalName = "D",
  data     = as.data.frame(raw)
)

# ---- 4) Harmonize, compute % and bounds ----
gsm <- gsm %>%
  mutate(
    Soil     = factor(Soil,     levels = soil_levels),
    Salinity = factor(Salinity, levels = salinity_levels),
    Moisture = factor(Moisture, levels = moisture_levels),
    Seeds    = as.numeric(Seeds),
    grs      = as.numeric(grs),
    mgt      = as.numeric(mgt)
  ) %>%
  filter(!is.na(Seeds), Seeds > 0, !is.na(grs), !is.na(mgt)) %>%
  mutate(
    grs       = pmin(grs, Seeds),
    germ_pct  = pmin(pmax(100 * (grs / Seeds), 0), 100)
  )

cat("\nQuick checks:\n")
cat("Rows:", nrow(gsm), "  Day cols detected:", paste(day_cols, collapse=", "), "\n")
print(table(Soil=gsm$Soil, Salinity=gsm$Salinity, Moisture=gsm$Moisture, useNA="ifany"))

# ---- 5) Basic GLMs ----
# A) Binomial GLM: germination proportion
mod_GP <- glm(
  cbind(grs, pmax(Seeds - grs, 0)) ~ Soil * Salinity * Moisture,
  data = gsm,
  family = binomial
)
cat("\n=== Binomial GLM for germination proportion ===\n"); print(summary(mod_GP))

# B) Gamma(log) GLM: mean germination time
mod_MGT <- glm(
  mgt ~ Soil * Salinity * Moisture,
  data = gsm,
  family = Gamma(link = "log")
)
cat("\n=== Gamma(log) GLM for MGT ===\n"); print(summary(mod_MGT))

# ---- 6) FIGURES: viridis boxplots (x=Moisture, fill=Salinity, facet=Soil) ----
time_labels <- c("Low"="Low","Medium"="Medium","High"="High")

fig_pct <- ggplot(gsm, aes(x = Moisture, y = germ_pct, fill = Salinity)) +
  geom_boxplot(outlier.alpha = 0.7, width = 0.7,
               position = position_dodge2(preserve = "single")) +
  scale_fill_viridis_d(name = "Salinity") +
  scale_y_continuous(limits = c(0,50), breaks = seq(0,50,10),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = time_labels, drop = FALSE) +
  labs(x = "Moisture Level", y = "Germination (%)") +
  facet_wrap(~ Soil, nrow = 1, labeller = labeller(Soil = soil_labels)) +  # <-- here
  theme_bw() +theme(legend.position = "right", 
                    axis.text.x = element_text(size = 11),
                    axis.title.x = element_text(size =12),
                    axis.title.y = element_text(size =12),
                    axis.text.y =  element_text(size =11),
                    legend.text = element_text(size = 11),
                    legend.title =  element_text(size =12))+ theme(legend.position = "bottom")



fig_pct




fig_mgt <- ggplot(gsm, aes(x = Moisture, y = mgt, fill = Salinity)) +
  geom_boxplot(outlier.alpha = 0.7, width = 0.7,
               position = position_dodge2(preserve = "single")) +
  scale_fill_viridis_d(name = "Salinity") +
  scale_x_discrete(labels = time_labels, drop = FALSE) +
  labs(x = "Moisture Level", y = "Days to germinate (MGT)") +
  facet_wrap(~ Soil, nrow = 1, labeller = labeller(Soil = soil_labels)) +  # <-- here
  theme_bw() + theme(legend.position = "bottom")


# Show & save
print(fig_pct); print(fig_mgt)
ggsave(file.path(out_dir, "LeeFigure3.tif"), fig_pct, dpi = 300)
ggsave(file.path(out_dir, "LeeSupplementalFigure7.tif"),fig_mgt, dpi = 300)

cat("\nDone. Figures saved to: ", normalizePath(out_dir), "\n")
