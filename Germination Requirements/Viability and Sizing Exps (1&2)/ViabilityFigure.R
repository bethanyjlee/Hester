rm(list = ls())

## Packages ##
need <- c("dplyr","ggplot2","scales","emmeans","multcompView","car")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales)
  library(emmeans); library(multcompView); library(car)})


## Set WD ##
setwd("~/R data/Chp1PWSeeds/ViabilityANDSizing")

## 1) Load + tidy to the actual sheet headers ##
dat <- read.csv("OctoberViability.csv", check.names = FALSE)

# Keep Early only - removed Later harvest from publication
dat <- dat %>%
  mutate(Harvest = as.character(Harvest)) %>%
  filter(!is.na(Harvest), grepl("^Early$", Harvest, ignore.case = TRUE))

# Merge Hester1/Hester2 → Hester
dat <- dat %>%
  mutate(`Seed Source` = ifelse(`Seed Source` %in% c("Hester1","Hester2"), "Hester", `Seed Source`))

# Basic checks/coercions
dat <- dat %>%
  mutate(Seeds = as.numeric(Seeds), Germ  = as.numeric(Germ))

# If Germ happened to be in % units, normalize (safe no-op if already 0–1)
if (max(dat$Germ, na.rm = TRUE) > 1) {
  dat <- dat %>% mutate(Germ = Germ / 100)}

# Build counts for modeling
dat <- dat %>%
  mutate( Germinated = pmax(round(Germ * Seeds), 0),
    Fail= pmax(Seeds - Germinated, 0))

## 2) Change site names to SiteIDs for Publication ##
site_map <- c(
  "MLWA"       = "1",
  "JettyRidge" = "2",
  "MoroCojo"   = "3",
  "SouthMarsh" = "4",  # Frequent
  "APC"        = "5",
  "Estrada"    = "6",
  "Kirby"      = "7",  # Frequent
  "Yampah"     = "8"   # Frequent
)

dat <- dat %>%
  mutate(
    SiteNum = dplyr::case_when(
      `Seed Source` %in% names(site_map) ~ unname(site_map[`Seed Source`]),
      `Seed Source` == "Hester"          ~ "Hester",
      TRUE                               ~ NA_character_),
    Inundation = dplyr::case_when(
      `Seed Source` %in% c("Kirby","Yampah","SouthMarsh") ~ "Frequent",
      TRUE                                                ~ "Infrequent")
  ) %>%
  filter(!is.na(SiteNum)) %>%
  mutate(
    Inundation = factor(Inundation, levels = c("Infrequent","Frequent")),
    # Order: Infrequent first (1,2,3,5,6,Hester) then Frequent (4,7,8)
    SiteNum    = factor(SiteNum, levels = c("1","2","3","5","6","Hester","4","7","8"))
  ) %>%
  droplevels()

## 3) Model (quasi-binomial) - Model Selection and checks done in Reference R file ##
glm_mod <- glm(
  cbind(Germinated, Fail) ~ SiteNum,
  data   = dat,
  family = quasibinomial("logit"))

emm <- emmeans(glm_mod, ~ SiteNum, type = "response")  # response-scale probs

## 4) Tukey letters (normalize contrast names to "A-B") ----
pairs_obj <- pairs(emm, adjust = "tukey")
pdat <- as.data.frame(summary(pairs_obj))

# Normalize contrast labels:
#  - turn "A / B" or "A/B" into "A-B"
#  - remove stray spaces
norm_contrast <- function(x) {
  x <- gsub("\\s*/\\s*", "-", x)   # slash to hyphen
  x <- gsub("\\s*-\\s*", "-", x)   # tidy spaced hyphens
  x <- gsub("\\s+", "", x)         # remove leftover spaces
  x
}
pdat$cn <- norm_contrast(pdat$contrast)

# Sanity check: every name must be exactly "lhs-rhs"
bad <- !grepl("^[^-]+-[^-]+$", pdat$cn)
if (any(bad)) {
  stop("Unexpected contrast labels after normalization: ",
       paste(unique(pdat$contrast[bad]), collapse = ", "))
}

# Named p-value vector for multcompView
pvec <- setNames(pdat$p.value, pdat$cn)

# Letters
lets <- multcompView::multcompLetters(pvec, threshold = 0.05)

letters_df <- data.frame(
  SiteNum = levels(dat$SiteNum),
  Letters = lets$Letters[levels(dat$SiteNum)],
  row.names = NULL)

## 5) Label placement summary for plotting ##
dt2 <- dat %>%
  group_by(SiteNum) %>%
  summarise(
    mean_germ = mean(Germ, na.rm = TRUE),
    max_germ  = max(Germ,  na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  left_join(letters_df, by = "SiteNum")

y_top <- min(1, max(dt2$max_germ, na.rm = TRUE) * 1.05 + 0.02)

## 6) Figure 1 ##
LeeFig1<-ggplot(dat, aes(x = SiteNum, y = Germ, fill = Inundation, cex.axis = 8)) +
  geom_boxplot(outlier.alpha = 0.5) +
  geom_text(
    data = dt2,
    inherit.aes = FALSE,
    aes(x = SiteNum, y = pmin(max_germ + 0.02, y_top), label = Letters),
    size = 5, vjust = 0) +
  labs(
    x = "Seed source",
    y = "Germination (%)",
    fill = "Inundation"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.55)) +
  scale_fill_manual(values = c("Infrequent" = "#1b9e77", "Frequent" = "#d95f02")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title = element_text(size =15),
        axis.text.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

ggsave("LeeFig1.tiff",LeeFig1, dpi=300)





data_summary <- dat %>%
  group_by(Inundation) %>%
 # filter(Date == max(Date)) %>%
  mutate(Germ) %>%
#  group_by(Seed.Source) %>%
  summarise(
    n         = n(),
    mean_germ = mean(Germ, na.rm = TRUE),
    sd_germ   = sd(Germ, na.rm = TRUE),
    se_germ   = sd_germ / sqrt(n)
  ) %>%
ungroup()

print(data_summary)
