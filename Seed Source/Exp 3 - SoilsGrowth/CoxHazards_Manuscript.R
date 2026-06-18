
rm(list=ls())
setwd("~/Hester/Seed Source/Exp 3 - SoilsGrowth")

library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)

#### 1) Read greenhouse sheet ####
greenhouse <- read.csv(
  "Elkhorn Greenhouse Data.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE)

#### Reshape to LONG format ####

# Step 1: Fix duplicate column names
names(greenhouse) <- c(
  "TagID", "Soil.Type", "TidalCat", "SiteID",
  "Date1", "Survivors", "Height1",
  "Date2", "Survivors2", "Height2",
  "Date3", "Survivors3", "Height3")

# Step 2: Reshape into long format
growth_long <- greenhouse %>%
  pivot_longer(cols = Date1:Height3,
               names_to = c(".value", "Timepoint"), 
               names_pattern = "(Date|Survivors|Height)([0-9]+)") %>%
  mutate(TagID = as.numeric(TagID),
         Height = as.numeric(Height),
         Survivors = as.numeric(Survivors)) %>%
  mutate(SoilType = factor(Soil.Type, levels = c("Hester", "Combo", "Potting Soil")),
         SiteID = factor(SiteID),
         Tidal = factor(TidalCat),
         Timepoint = factor(Timepoint,
                            levels = c("1","2","3"),
                            labels = c("6/12", "7/22", "8/21")),
         Date = factor(Date))

#### 2) Read seed source sheet ####
seedsource <- read.csv("FullGermElkhornSeedSource.csv", stringsAsFactors = FALSE)

## check names if needed
# names(seedsource)

## make sure TagID exists and matches type
dat <- greenhouse %>%
  mutate(
    TagID = as.numeric(TagID),
    SoilType = as.character(Soil.Type),
    `Survivors` = as.numeric(`Survivors`),
    `Survivors2` = as.numeric(`Survivors2`),
    `Survivors3` = as.numeric(`Survivors3`))

seedsource <- seedsource %>%
  mutate(TagID = as.numeric(TagID))

#### 3) Join seed source #### 
dat2 <- dat %>%
  left_join(seedsource, by = "TagID")



#### 4) Clean survival counts #### If start is missing, assume 5 seedlings
dat2 <- dat2 %>%
  mutate(
    Start = Survivors,
    S1 = Start,
    S2 = `Survivors2`,
    S3 = `Survivors3`,
    S4 = `Survivors3`
  )

## necessary safety checks
dat2 <- dat2 %>%
  mutate(
    S1 = pmax(S1, 0),
    S2 = pmax(S2, 0),
    S3 = pmax(S3, 0),
    S4 = pmax(S4, 0)
  )

## enforce non-increasing survival if there are data entry hiccups
dat2 <- dat2 %>%
  rowwise() %>%
  mutate(
    S2 = min(S1, S2, na.rm = TRUE),
    S3 = min(S2, S3, na.rm = TRUE),
    S4 = min(S3, S4, na.rm = TRUE)
  ) %>%
  ungroup()

#### 5) Convert census counts into individual seedling survival times ####
## times are days since first count date (6/12)
time_0 <- 0
time_1 <- as.numeric(as.Date("2024-07-22") - as.Date("2024-06-12"))
time_2 <- as.numeric(as.Date("2024-08-21") - as.Date("2024-06-12"))
time_3 <- as.numeric(as.Date("2024-09-24") - as.Date("2024-06-12"))


#### 6) Build individual survival dataset #### 
# the big step - patience is needed #

surv_data <- dat2 %>%
  rowwise() %>%
  do({
    row <- .
    
    # pull values safely
    s1 <- as.numeric(row$S1)
    s2 <- as.numeric(row$S2)
    s3 <- as.numeric(row$S3)
    s4 <- as.numeric(row$S4)
    
    # replace NA with 0
    s1 <- ifelse(is.na(s1), 0, s1)
    s2 <- ifelse(is.na(s2), 0, s2)
    s3 <- ifelse(is.na(s3), 0, s3)
    s4 <- ifelse(is.na(s4), 0, s4)
    
    # enforce non-increasing again (extra safety)
    s2 <- min(s1, s2)
    s3 <- min(s2, s3)
    s4 <- min(s3, s4)
    
    # calculate deaths
    d1 <- max(0, s1 - s2)
    d2 <- max(0, s2 - s3)
    d3 <- max(0, s3 - s4)
    
    # ensure integers
    d1 <- as.integer(d1)
    d2 <- as.integer(d2)
    d3 <- as.integer(d3)
    s4 <- as.integer(s4)
    
    # build dataframe safely
    df <- data.frame(
      time = c(
        rep(time_1, d1),
        rep(time_2, d2),
        rep(time_3, d3),
        rep(time_3, s4)
      ),
      status = c(
        rep(1, d1),
        rep(1, d2),
        rep(1, d3),
        rep(0, s4)
      )
    )
    
    # only add metadata if rows exist
    if (nrow(df) > 0) {
      df$SiteID <- row$SiteID.x
      df$Tidal <- row$TidalCat.x
      df$SoilType <- row$Soil.Type.x
      df$TagID <- row$TagID
    }
    
    df
  }) %>%
  ungroup()


#### Run Baseline Models ####

# Site-level model
cox_site <- survfit(Surv(time, status) ~ SiteID, data = surv_data)
summary(cox_site)

# Tidal-level model
tidal_site <- survfit(Surv(time, status) ~ Tidal, data = surv_data)
summary(tidal_site)

#Soil model
cox_soil <- coxph(Surv(time, status) ~ SoilType, data = surv_data)
summary(cox_soil)

#Interaction model - Site*Soil
cox_soil_site_int <- coxph(Surv(time, status) ~ SoilType * SiteID, data = surv_data)
summary(cox_soil_site_int)

# Model - Site n Soil
cox_soil_site <- coxph(Surv(time, status) ~ SoilType + SiteID, data = surv_data)
summary(cox_soil_site)

#Interaction model - Tidal*Soil
cox_soil_tide_int <- coxph(Surv(time, status) ~ SoilType * Tidal, data = surv_data)
summary(cox_soil_tide_int)

# Model - Tidal n Soil
cox_soil_tide <- coxph(Surv(time, status) ~ SoilType + Tidal, data = surv_data)
summary(cox_soil_tide)

# Make SiteID Random
Cox_soil_randSite<-coxph(Surv(time, status) ~ SoilType + frailty(SiteID), data =surv_data)
summary(Cox_soil_randSite)


#### check assumptions ####
ph_test <- cox.zph(Cox_soil_randSite)
print(ph_test)
plot(ph_test)


fit_soil <- survfit(cox_soil, newdata = data.frame(
  SoilType = factor(c("Hester", "Combo", "Potting Soil"),
                    levels = c("Hester", "Combo", "Potting Soil"))))

ggsurvplot(
  fit_soil,
  data = surv_data,
  conf.int = TRUE,
  legend.title = "Soil Type",
  legend.labs = c("Hester", "Combo", "Potting Soil"),
  xlab = "Days since 6/12",
  ylab = "Predicted survival probability",
  ggtheme = theme_bw(),
  risk.table = TRUE)


#ggforest(cox_full, data = surv_dat)

km_soil <- survfit(Surv(time, status) ~ SoilType, data = surv_data)

ggsurvplot(
  km_soil,
  data = surv_data,
  conf.int = TRUE,
  pval = TRUE,
  #risk.table = TRUE,
  xlab = "Days since 6/12",
  ylab = "Survival probability",
  legend.title = "Soil Type",
  ggtheme = theme_bw())

ggsave(filename = "Figure6.tif", dpi = 300, path = "Figures")


surv_data$SeedSource <- factor(surv_data$SiteID)

km_seed <- survfit(Surv(time, status) ~ SiteID, data = surv_data)

ggsurvplot(
  km_seed,
  data = surv_data,
  conf.int = FALSE,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Days since 6/12",
  ylab = "Survival probability",
  legend.title = "Seed Source",
  ggtheme = theme_bw())

ggsurvplot(
  km_seed,
  data = surv_data,
  conf.int = FALSE,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Days since 6/12",
  ylab = "Survival probability",
  legend.title = "Seed Source",
  ggtheme = theme_bw(),
  legend = "bottom")

ggsurvplot(
  km_seed,
  data = surv_data,
  conf.int = FALSE,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Days since 6/12",
  ylab = "Survival probability",
  legend.title = "Seed Source",
  ggtheme = theme_bw(),
  legend = "bottom")



# Define custom soil colors
soil_colors <- c(
  "Hester" = "#999333",
  "Combo" = "#aa4499",
  "Potting Soil" = "#661100")

# Replot model-based survival (fit_soil) with colors
ggsurvplot(
  fit_soil,
  data = surv_data,
  conf.int = TRUE,
 # risk.table = TRUE,
  palette = soil_colors,   # << key line
  legend.title = "Soil Type",
  legend.labs = c("Hester", "Combo", "Potting Soil"),
  xlab = "Days since 6/12",
  ylab = "Predicted survival probability",
  ggtheme = theme_bw())


#### Kaplan-Meier by Tidal ####
km_tidal <- survfit(Surv(time, status) ~ Tidal, data = surv_data)

ggsurvplot(
  km_tidal,
  data = surv_data,
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  xlab = "Days since 6/12",
  ylab = "Survival probability",
  legend.title = "Tidal Category",
  ggtheme = theme_bw())

#### Cox model by Tidal ####
cox_tidal_only <- coxph(Surv(time, status) ~ Tidal, data = surv_data)
summary(cox_tidal_only)

#### Optional: proportional hazards check ####
ph_tidal <- cox.zph(cox_tidal_only)
print(ph_tidal)
plot(ph_tidal)


### Relevel SiteID so your preferred reference is used
# replace "Restored1" with whichever site you want as the reference
surv_data$SiteID <- relevel(
  factor(surv_data$SiteID),
  ref = "Restored1")

### Cox model with SiteID only
cox_site_hr <- coxph(
  Surv(time, status) ~ SiteID,
  data = surv_data)

summary(cox_site_hr)

### Extract hazard ratios + confidence intervals
site_hr <- data.frame(
  Site = rownames(summary(cox_site_hr)$coefficients),
  HR = summary(cox_site_hr)$conf.int[, "exp(coef)"],
  LowerCI = summary(cox_site_hr)$conf.int[, "lower .95"],
  UpperCI = summary(cox_site_hr)$conf.int[, "upper .95"],
  pvalue = summary(cox_site_hr)$coefficients[, "Pr(>|z|)"])

### Clean site names
site_hr$Site <- gsub("SiteID", "", site_hr$Site)

### Manually set site order
site_order <- c(
  "Muted1",
  "Muted2",
  "Muted3",
  "Muted4",
  "Muted5",
  "Restored1",
  "Tidal2",
  "Tidal4",
  "Tidal5")

### Apply ordering
site_hr$Site <- factor(
  site_hr$Site,
  levels = rev(site_order)))

### Add reference site back into dataframe
ref_row <- data.frame(
  Site = "Restored1",
  HR = 1,
  LowerCI = 1,
  UpperCI = 1,
  pvalue = NA)

site_hr <- rbind(site_hr, ref_row)


### Hazard ratio plot
ggplot(site_hr, aes(x = Site, y = HR)) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "gray40") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerCI,
                    ymax = UpperCI),
                width = 0.2,
                linewidth = 0.8) +
  coord_flip() +
  theme_bw() +
  labs(x = "Seed Source Site",
    y = "Hazard Ratio (95% CI)",  ) +
  theme(panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(), text = element_text(size = 13))


### Summary table 

#### Final Survival Summary Table by Seed Source ####

survival_summary <- dat2 %>%
  group_by(SiteID.x) %>%
  summarise(
    Starting_Seedlings = sum(S1, na.rm = TRUE),
    Final_Survivors = sum(S4, na.rm = TRUE),
    Survival_Proportion = Final_Survivors / Starting_Seedlings
  ) %>%
  arrange(SiteID.x)

print(survival_summary)


anova(cox_soil_tide, cox_soil_tide_int, test = "LRT")

cox_soil_only <- coxph(Surv(time, status) ~ SoilType, data = surv_data)

anova(cox_soil_only, cox_soil_tide, test = "LRT")

drop1(Cox_soil_randSite, test = "Chisq")



####emmeans ####

library(emmeans)

emm <- emmeans(
  cox_soil_tide_int,
  ~ SoilType | Tidal
)

pairs(emm)
