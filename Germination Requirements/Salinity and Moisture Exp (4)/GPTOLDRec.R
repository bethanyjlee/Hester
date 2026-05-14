fulldata<-read.csv('MegaExp.csv')

long<- fulldata %>%
  pivot_longer(
    cols = starts_with("D"),
    names_to= "Date")

long<- long %>%
  mutate(Date= as.numeric(gsub("D","", Date)))

long<-long %>%
  mutate(value = as.double(value))

long<- long %>%
  mutate(Salinity = as.character(Salinity))

dataready1<- long %>%
  mutate(across(c(Soil,Salinity,Moisture),as.factor))

plottiming <-dataready1 %>%
  group_by(Date, Soil,Salinity,Moisture) %>%
  summarise(germs = value/ Seeds)

str(plottiming)

summ<- plottiming %>%
  group_by(Date,Soil,Salinity,Moisture) %>%
  summarize(mean = mean(germs))

summ$size_f = factor(summ$Moisture, levels = c('Low','Medium','High'))


query1<-summ$Soil=='PottingSoil'
index1<-which(query1)
Potting<-summ[index1,]
Mimic<-summ[-index1,]


summ %>%
  ggplot(aes(x=Date, y=mean)) +
  geom_point() +
  stat_smooth(method = "lm")+
  labs(x='Days After Sown', y= "Germination of Salicornia",
       title = "Response of Salicornia Germination") +
  scale_shape_manual(values=c(0,15))+
  scale_color_manual(values=c('red','blue'))+
  stat_cor(label.y = 0.4) +
  stat_regline_equation(label.y = 0.3)+
  facet_grid(size_f~ Soil + Salinity)+
  theme_bw()

?stat_cor
#### try to determine interactions ####

model<-aov(germs~Soil + Salinity + Moisture, data = plottiming)
summary(model)
plot(model) ## fits decently?

## try a vif first
vif(lm(germs~Soil + Salinity + Moisture, data = plottiming))

#Salinity and Moisture
interaction.plot(plottiming$Salinity, plottiming$Moisture, plottiming$germs)
#lines cross = interaction
#Soil and Moisture
interaction.plot(plottiming$Soil, plottiming$Moisture, plottiming$germs)
# lines cross = interaction
interaction.plot(plottiming$Soil, plottiming$Salinity, plottiming$germs)


modelint<-lm(germs~Soil*Moisture*Salinity, data = plottiming)
summary(modelint)
plot(modelint)


#### germination and emergence GerminaR ####
#### new datasheet to show only new growth per day ####

germsum<-read.csv('Mega2New.csv')

query1<-germsum$Soil=='PottingSoil'
index1<-which(query1)
Potting<-germsum[index1,]
Mimic<-germsum[-index1,]

gsm<-ger_summary(SeedN = "Seeds", evalName = "D", data=germsum)
summary(gsm)


summary(gsm)

hester<-ger_summary(SeedN = "Seeds", evalName = "D", data=Mimic)
potting<-ger_summary(SeedN = "Seeds", evalName = "D", data=Potting)

query2<-germsum$Salinity=='Saltwater'
index2<-which(query2)
Fresh<-germsum[index2,]
Salt<-germsum[-index2,]

fresh<-ger_summary(SeedN = "Seeds", evalName = "D", data=Fresh)
salt<-ger_summary(SeedN = "Seeds", evalName = "D", data=Salt)

query3<-germsum$Moisture=='Low'
index3<-which(query3)
LowM<-germsum[index3,]

query4<-germsum$Moisture=='Medium'
index4<-which(query4)
MediumM<-germsum[index4,]

query5<-germsum$Moisture=='High'
index5<-which(query5)
HighM<-germsum[index5,]


low<-ger_summary(SeedN = "Seeds", evalName = "D", data=LowM)
med<-ger_summary(SeedN = "Seeds", evalName = "D", data=MediumM)
high<-ger_summary(SeedN = "Seeds", evalName = "D", data=HighM)

summary(hester)
summary(potting)
summary(salt)
summary(fresh)
summary(low)
summary(med)
summary(high)

#### ANOVA for mgt ####

av <- aov(mgt ~ Salinity*Moisture, data = gsm)
summary(av)

HMgsm<- aov(mgt~Salinity*Moisture, data = hester)
summary(HMgsm)
PSgsm<-aov(mgt~Salinity*Moisture, data = potting)
summary(PSgsm)

# mean comparison test

mc_mgt <- ger_testcomp(aov = av, comp = c("Salinity", "Moisture"))


mc_mgt$table %>% 
  kable(caption = "Mean germination time comparison")


mgt <- mc_mgt$table %>% 
  fplot(data = .
        , type = "bar" 
        , x = "Salinity"
        , y = "mgt"
        , group = "Moisture"
        , ylimits = c(0,15, 1)
        , ylab = "Mean germination time (days)"
        , xlab = "Salinity"
        , glab = "Moisture Levels"
        , sig = "sig"
        , error = "ste"
        , color = T)

mgt


## potting soil only!! ##

# mean comparison test

PSmc_mgt <- ger_testcomp(aov = PSgsm, comp = c("Salinity", "Moisture"))

# data result

PSmc_grs


mc_mgt$table %>% 
  kable(caption = "Mean germination time comparison")


PSmgt <- PSmc_mgt$table %>% 
  fplot(data = .
        , type = "bar" 
        , x = "Salinity"
        , y = "mgt"
        , group = "Moisture"
        , ylimits = c(0,15, 1)
        , ylab = "Mean germination time (days)"
        , xlab = "Salinity"
        , glab = "Moisture Levels"
        , sig = "sig"
        , color = T)

PSmgt

?fplot

### hester soil time


# mean comparison test

HMmc_mgt <- ger_testcomp(aov = HMgsm, comp = c("Salinity", "Moisture"))

# data result


mc_mgt$table %>% 
  kable(caption = "Mean germination time comparison")


HMmgt <- HMmc_mgt$table %>% 
  fplot(data = .
        , type = "bar" 
        , x = "Salinity"
        , y = "mgt"
        , group = "Moisture"
        , ylimits = c(0,15, 1)
        , ylab = "Mean germination time (days)"
        , xlab = "Salinity"
        , glab = "Moisture Levels"
        , sig = "sig"
        , error = "ste"
        , color = T)

HMmgt


#### ANOVA for grs ####

## all soils

avg <- aov(grs ~ Salinity*Moisture, data = gsm)
summary(avg)

HMavg<- aov(grs~Salinity*Moisture, data = hester)
summary(HMavg)
PSavg<-aov(grs~Salinity*Moisture, data = potting)
summary(PSavg)

HMavg$Model$grs$MoistureF = factor(HMavg$Model$grs$Moisture, levels = c('Low','Medium','High'))
PSavg$Model$grs$MoistureF = factor(PSavg$Model$grs$Moisture, levels = c('Low','Medium','High'))

# mean comparison test

mc_grs <- ger_testcomp(aov = avg, comp = c("Salinity", "Moisture"))

mc_grs$table$MoistureF = factor(mc_grs$table$Moisture, levels = c('Low','Medium','High'))

grc <- mc_grs$table %>% 
  fplot(data = .
        , type = "bar" 
        , x = "Salinity"
        , y = "grs"
        , group = "MoistureF"
        , ylab = "Mean Germinated Seeds"
        , xlab = "Salinity"
        , glab = "Moisture Levels"
        , sig = "sig"
        , error = "ste"
        , color = T)

grc


## potting soil only!! ##

summary(PSavg)

# mean comparison test

PSmc_grs <- ger_testcomp(aov = PSavg, comp = c("Salinity", "Moisture"))

PSmc_grs$table$MoistureF = factor(mc_grs$table$Moisture, levels = c('Low','Medium','High'))

PSgrc <- PSmc_grs$table %>% 
  fplot(data = .
        , type = "bar" 
        , x = "Salinity"
        , y = "grs"
        , group = "MoistureF"
        , ylab = "Mean Germinated Seeds"
        , xlab = "Salinity"
        , glab = "Potting Soil                   Moisture Levels"
        , sig = "sig"
        , error = "ste"
        , color = T)



PSgrc

?fplot

### hester soil time

summary(HMavg)

# mean comparison test

HMmc_grs <- ger_testcomp(aov = HMavg, comp = c("Salinity", "Moisture"))

HMmc_grs$table$MoistureF = factor(mc_grs$table$Moisture, levels = c('Low','Medium','High'))

# data result


mc_mgt$table %>% 
  kable(caption = "Mean germination time comparison")


HMgrs <- HMmc_grs$table %>% 
  fplot(data = .
        , type = "bar" 
        , x = "Salinity"
        , y = "grs"
        , group = "MoistureF"
        , ylab = "Mean Germinated Seeds"
        , xlab = "Salinity"
        , glab = "Hester Mimic                   Moisture Levels"
        , sig = "sig"
        , error = "ste"
        , color = T)

HMgrs




#### Treatment separate ANOVA for Germ ####

PSsaltlow <- Potting %>%
  filter(Salinity == "30") %>%
  filter(Moisture == "Low")

PSsaltmed <- Potting %>%
  filter(Salinity == "30") %>%
  filter(Moisture == "Medium")

PSsalthigh <- Potting %>%
  filter(Salinity == "30") %>%
  filter(Moisture == "High")


PSfreshlow <- Potting %>%
  filter(Salinity == "0") %>%
  filter(Moisture == "Low")

PSfreshmed <- Potting %>%
  filter(Salinity == "0") %>%
  filter(Moisture == "Medium")

PSfreshhigh <- Potting %>%
  filter(Salinity == "0") %>%
  filter(Moisture == "High")





HMsaltlow <- Mimic %>%
  filter(Salinity == "30") %>%
  filter(Moisture == "Low")

HMsaltmed <- Mimic %>%
  filter(Salinity == "30") %>%
  filter(Moisture == "Medium")

HMsalthigh <- Mimic %>%
  filter(Salinity == "30") %>%
  filter(Moisture == "High")


HMfreshlow <- Mimic %>%
  filter(Salinity == "0") %>%
  filter(Moisture == "Low")

HMfreshmed <- Mimic %>%
  filter(Salinity == "0") %>%
  filter(Moisture == "Medium")

HMfreshhigh <- Mimic %>%
  filter(Salinity == "0") %>%
  filter(Moisture == "High")



PSSLow<-ger_summary(SeedN = "Seeds", evalName = "D", data=PSsaltlow)
PSSMed<-ger_summary(SeedN = "Seeds", evalName = "D", data=PSsaltmed)
PSSHigh<-ger_summary(SeedN = "Seeds", evalName = "D", data=PSsalthigh)

summary(PSSLow)
summary(PSSMed)
summary(PSSHigh)


PSFLow<-ger_summary(SeedN = "Seeds", evalName = "D", data=PSfreshlow)
PSFMed<-ger_summary(SeedN = "Seeds", evalName = "D", data=PSfreshmed)
PSFHigh<-ger_summary(SeedN = "Seeds", evalName = "D", data=PSfreshhigh)

summary(PSFLow)
summary(PSFMed)
summary(PSFHigh)




HMSLow<-ger_summary(SeedN = "Seeds", evalName = "D", data=HMsaltlow)
HMSMed<-ger_summary(SeedN = "Seeds", evalName = "D", data=HMsaltmed)
HMSHigh<-ger_summary(SeedN = "Seeds", evalName = "D", data=HMsalthigh)

summary(HMSLow)
summary(HMSMed)
summary(HMSHigh)


HMFLow<-ger_summary(SeedN = "Seeds", evalName = "D", data=HMfreshlow)
HMFMed<-ger_summary(SeedN = "Seeds", evalName = "D", data=HMfreshmed)
HMFHigh<-ger_summary(SeedN = "Seeds", evalName = "D", data=HMfreshhigh)

summary(HMFLow)
summary(HMFMed)
summary(HMFHigh)


#### plots relevant for MS ####

#### Potting soil

# mean comparison test

PS_grs <- ger_testcomp(aov = PSavg, comp = c("Moisture","Salinity"))


PS_grs$table$SalinityF = factor(PS_grs$table$Salinity, levels = c('Low','Medium','High'))
PS_grs$table$MoistureF = factor(PS_grs$table$Moisture, levels = c('Freshwater','Saltwater'))

table<-PS_grs$table


ggplot(data = table, aes(x=Salinity, y=grs, fill = Moisture)) +
  geom_col(position = position_dodge2(preserve = 'single')) + 
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Seeds Germinated in High Nutrient, Potting Soil", x = "Addition to Soils",
       y = "Number of Germinated Seeds", fill = "Moisture Levels")



### hester mimic soil 

HM_grs <- ger_testcomp(aov = HMavg, comp = c("Moisture","Salinity"))


HM_grs$table$SalinityF = factor(HM_grs$table$Salinity, levels = c('Low','Medium','High'))
HM_grs$table$MoistureF = factor(HM_grs$table$Moisture, levels = c('Freshwater','Saltwater'))

table2<-HM_grs$table


ggplot(data = table2, aes(x=Salinity, y=grs, fill = Moisture)) +
  geom_col(position = position_dodge2(preserve = 'single')) + 
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Seeds Germinated in Low Nutrient, Clay Soil", x = "Addition to Soils",
       y = "Number of Germinated Seeds", fill = "Moisture Levels")


#### updated ####

potting$MoistureF = factor(potting$Moisture, levels = c('Low','Medium','High'))
potting$SalinityF = factor(potting$Salinity, levels = c('Freshwater','Saltwater'))

PSGRS<-ggplot(data = potting, aes(x = Salinity, y = grs, fill = MoistureF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) +  # Bar width adjustment
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Seeds Germinated in High Nutrient, Potting Soil", 
       x = "Moisture Addition to Soils", 
       y = "Number of Germinated Seeds", 
       fill = "Moisture Levels") +
  scale_y_continuous(limits = c(0,35))+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw()


PSMGT<-ggplot(data = potting, aes(x = Salinity, y = mgt, fill = MoistureF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) +  # Bar width adjustment
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Days to Germinate in High Nutrient, Potting Soil", 
       x = "Moisture Addition to Soils", 
       y = "Days to Germinate", 
       fill = "Moisture Levels") +
  scale_y_continuous(limits = c(0,16))+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw()



hester$MoistureF = factor(hester$Moisture, levels = c('Low','Medium','High'))
hester$SalinityF = factor(hester$Salinity, levels = c('Freshwater','Saltwater'))

HMGRS<-ggplot(data = hester, aes(x = Salinity, y = grs, fill = MoistureF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) +  # Bar width adjustment
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Seeds Germinated in Nutrient Poor, Hester Soil", 
       x = "Moisture Addition to Soils", 
       y = "Number of Germinated Seeds", 
       fill = "Moisture Levels") +
  scale_y_continuous(limits = c(0,35))+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw()


HMMGT<-ggplot(data = hester, aes(x = Salinity, y = mgt, fill = MoistureF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) +  # Bar width adjustment
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Days to Germinate in Nutrient Poor, Hester Soil", 
       x = "Moisture Addition to Soils", 
       y = "Days to Germinate", 
       fill = "Moisture Levels") +
  scale_y_continuous(limits = c(0,16))+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw()


library(gridExtra)

grid.arrange(PSGRS,HMGRS,PSMGT,HMMGT,ncol=2,nrow=2)
plot(HMGRS)
plot(HMMGT)
plot(PSGRS)
plot(PSMGT)

#### 8/11/25 ####
# PACKAGES ----
library(tidyverse)
library(car)         # Type-II/III tests
library(emmeans)     # marginal means & contrasts
library(performance) # overdispersion checks
library(DHARMa)      # residual diagnostics
library(broom)       # tidy outputs
library(GerminaR)    # mgt & grs summaries
library(ggpubr)      # optional, for arranging plots



# 0) HELPER: set factor orders you want in figures/tables
moisture_levels <- c("Low","Medium","High")
salinity_levels <- c("Freshwater","Saltwater")   # 0 ppt vs 30 ppt
soil_levels     <- c("PottingSoil","HesterMimic")# rename if needed

# 1) LOAD & SHAPE DAILY CUMULATIVE DATA (FINAL GERMINATION) ----
raw <- read.csv("MegaExp.csv", na.strings = c("", "NA"))

# Identify day columns (cumulative counts)
day_cols <- grep("^D\\d+$", names(raw), value = TRUE)
stopifnot(length(day_cols) > 0)

# Ensure expected columns exist
needed <- c("Soil","Salinity","Moisture","Seeds")
missing_needed <- setdiff(needed, names(raw))
if (length(missing_needed)) {
  stop(paste("Missing columns in MegaExp.csv:", paste(missing_needed, collapse=", ")))
}

# Optional: bring in a replicate id if present (Tray/Rep/etc.)
rep_col <- intersect(names(raw), c("Tray","Rep","Replicate","UnitID"))
rep_col <- if (length(rep_col)) rep_col[1] else NA_character_

# Long format for cumulative daily counts
long_cum <- raw %>%
  pivot_longer(all_of(day_cols), names_to = "Day", values_to = "cum_germ") %>%
  mutate(
    Day = as.numeric(gsub("^D","", Day)),
    cum_germ = as.numeric(cum_germ)
  )

# Build a unique unit id
if (is.na(rep_col)) {
  long_cum <- long_cum %>%
    group_by(Soil, Salinity, Moisture) %>%
    mutate(Unit = cur_group_id()) %>%
    ungroup()
} else {
  long_cum <- long_cum %>% rename(Unit = !!rep_col)
}

# Keep the final cumulative germination per unit = max across days
final <- long_cum %>%
  group_by(Unit, Soil, Salinity, Moisture) %>%
  summarise(
    Seeds = suppressWarnings(max(raw$Seeds[match(Unit, if (is.na(rep_col)) {
      # if we constructed Unit, retrieve seeds within grouping
      which(raw$Soil==first(Soil) & raw$Salinity==first(Salinity) & raw$Moisture==first(Moisture))[1]
    } else {
      match(Unit, raw[[rep_col]])
    })], na.rm = TRUE)),
    final_germ = max(cum_germ, na.rm = TRUE),
    .groups = "drop"
  )

# If the previous Seeds retrieval feels brittle (depends on your file layout), prefer:
# -> Pull Seeds directly from raw in the long join instead of matching:
# final <- long_cum %>%
#   left_join(raw %>% select(all_of(c(rep_col, "Soil","Salinity","Moisture","Seeds"))),
#             by = if (is.na(rep_col)) c("Soil","Salinity","Moisture") else c("Unit" = rep_col)) %>%
#   group_by(Unit, Soil, Salinity, Moisture, Seeds) %>%
#   summarise(final_germ = max(cum_germ, na.rm = TRUE), .groups = "drop")

# Clean and constrain
final_clean <- final %>%
  mutate(
    Seeds = as.integer(Seeds),
    final_germ = pmax(pmin(as.integer(round(final_germ)), Seeds), 0),
    Soil     = factor(Soil,     levels = soil_levels),
    Salinity = factor(Salinity, levels = salinity_levels),
    Moisture = factor(Moisture, levels = moisture_levels)
  ) %>%
  drop_na(Soil, Salinity, Moisture, Seeds, final_germ) %>%
  filter(Seeds > 0)



# 2) BINOMIAL GLM FOR FINAL GERMINATION PROPORTION ----
# Full factorial: Soil * Salinity * Moisture
m_binom <- glm(
  cbind(final_germ, Seeds - final_germ) ~ Soil * Salinity * Moisture,
  family = binomial,
  data   = final_clean
)

# Overdispersion check and fallback to quasibinomial if needed
od <- performance::check_overdispersion(m_binom)
print(od)

m_final <- if (od$dispersion_ratio > 1.2 && od$p_value < 0.05) {
  message("Refitting with quasibinomial due to overdispersion...")
  glm(
    cbind(final_germ, Seeds - final_germ) ~ Soil * Salinity * Moisture,
    family = quasibinomial,
    data   = final_clean
  )
} else m_binom

# Type-II tests (change to type="III" if you set sum-to-zero contrasts)
Anova(m_final, type = 2)

# DHARMa residual diagnostics (works with binomial; with quasibinomial it’s still informative)
sim <- simulateResiduals(m_final, plot = FALSE)
plot(sim)  # inspect QQ, residuals vs. fitted, etc.

# 3) EMMEANS, PAIRWISE CONTRASTS, PREDICTIONS ----
emm <- emmeans(m_final, ~ Soil * Salinity * Moisture, type = "response")
# Pairwise by Moisture within Soil×Salinity (edit as you like)
pairs(emm, by = c("Soil","Salinity"), adjust = "tukey")

# Predicted probabilities and 95% CIs for plotting
pred_df <- as.data.frame(emm)

# Plot predicted germination (%) with CIs
p_pred <- ggplot(pred_df,
                 aes(Moisture, prob, ymin = asymp.LCL, ymax = asymp.UCL,
                     color = Salinity, group = Salinity)) +
  geom_point(position = position_dodge(width = 0.35)) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.35)) +
  geom_line(position = position_dodge(width = 0.35)) +
  facet_wrap(~ Soil, nrow = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Moisture level",
       y = "Predicted germination (%)",
       color = "Salinity") +
  theme_bw()
print(p_pred)

# 4) GERMINATION TIME (mgt) & SPEED (grs) WITH GerminaR ----
germsum <- read.csv("Mega2New.csv", na.strings = c("", "NA"))




#### 9/11/25 ####

# PACKAGES ----
library(tidyverse)
library(car)         # Type-II/III tests
library(emmeans)     # marginal means & contrasts
library(performance) # overdispersion checks
library(DHARMa)      # residual diagnostics
library(broom)       # tidy outputs
library(GerminaR)    # mgt & grs summaries (optional for reporting)
library(ggpubr)      # optional, for arranging plots

# FACTOR ORDERS ----
moisture_levels <- c("Low","Medium","High")
salinity_levels <- c("Freshwater","Saltwater")
soil_levels     <- c("PottingSoil","HesterMimic")

# 1) LOAD & SHAPE DAILY CUMULATIVE DATA (FINAL GERMINATION) ----
raw <- read.csv("MegaExp.csv", na.strings = c("", "NA"))

day_cols <- grep("^D\\d+$", names(raw), value = TRUE)
stopifnot(length(day_cols) > 0)

needed <- c("Soil","Salinity","Moisture","Seeds")
missing_needed <- setdiff(needed, names(raw))
if (length(missing_needed)) {
  stop(paste("Missing columns in MegaExp.csv:", paste(missing_needed, collapse=", ")))
}

# Prefer a clean long join that carries Seeds forward reliably
rep_col <- intersect(names(raw), c("Tray","Rep","Replicate","UnitID"))
rep_col <- if (length(rep_col)) rep_col[1] else NA_character_

long_cum <- raw %>%
  pivot_longer(all_of(day_cols), names_to = "Day", values_to = "cum_germ") %>%
  mutate(
    Day = as.numeric(gsub("^D","", Day)),
    cum_germ = as.numeric(cum_germ)
  )

if (is.na(rep_col)) {
  # No explicit replicate column: carry Seeds and create a Unit id
  long_cum <- raw %>%
    select(Soil, Salinity, Moisture, Seeds, all_of(day_cols)) %>%
    pivot_longer(all_of(day_cols), names_to = "Day", values_to = "cum_germ") %>%
    mutate(
      Day = as.numeric(sub("^D","", Day)),
      cum_germ = as.numeric(cum_germ)
    ) %>%
    group_by(Soil, Salinity, Moisture, Seeds) %>%
    mutate(Unit = cur_group_id()) %>%
    ungroup()
} else {
  # Replicate column present: rename to Unit and carry Seeds
  long_cum <- raw %>%
    select(all_of(rep_col), Soil, Salinity, Moisture, Seeds, all_of(day_cols)) %>%
    rename(Unit = !!rep_col) %>%
    pivot_longer(all_of(day_cols), names_to = "Day", values_to = "cum_germ") %>%
    mutate(
      Day = as.numeric(sub("^D","", Day)),
      cum_germ = as.numeric(cum_germ)
    )
}

# --- FINAL CUMULATIVE (per Unit) ---
final_clean <- long_cum %>%
  group_by(Unit, Soil, Salinity, Moisture, Seeds) %>%
  summarise(final_germ = max(cum_germ, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Seeds      = as.integer(Seeds),
    final_germ = pmax(pmin(as.integer(round(final_germ)), Seeds), 0),
    Soil       = factor(Soil,     levels = soil_levels),
    Salinity   = factor(Salinity, levels = salinity_levels),
    Moisture   = factor(Moisture, levels = moisture_levels)
  ) %>%
  drop_na(Soil, Salinity, Moisture, Seeds, final_germ) %>%
  filter(Seeds > 0)

# 2) GLM-ONLY WORKFLOW WITH OVERDISPERSION CHECK ----
# Step 1: Fit binomial GLM
m_binom <- glm(
  cbind(final_germ, Seeds - final_germ) ~ Soil * Salinity * Moisture,
  family = binomial,
  data   = final_clean
)

# Step 2: Check overdispersion on the binomial fit
od <- performance::check_overdispersion(m_binom)
print(od)

# Step 3: If overdispersed, refit with quasibinomial; otherwise keep binomial
m_final <- if (isTRUE(od$dispersion_ratio > 1.2) && isTRUE(od$p_value < 0.05)) {
  message("Refitting with quasibinomial due to overdispersion...")
  glm(
    cbind(final_germ, Seeds - final_germ) ~ Soil * Salinity * Moisture,
    family = quasibinomial,
    data   = final_clean
  )
} else {
  m_binom
}

# 3) INFERENCE, DIAGNOSTICS, EMMs, PLOTS ----
# Type-II tests (use type = "III" if you set sum-to-zero contrasts)
Anova(m_final, type = 2)

# DHARMa diagnostics
sim <- simulateResiduals(m_final, plot = FALSE)
plot(sim)

# EMMs on response (probability) scale and Tukey-adjusted pairwise contrasts
emm <- emmeans(m_final, ~ Soil * Salinity * Moisture, type = "response")
pairs(emm, by = c("Soil","Salinity"), adjust = "tukey")

# Predicted probabilities for plotting
pred_df <- as.data.frame(emm)

p_pred <- ggplot(pred_df,
                 aes(Moisture, prob, ymin = asymp.LCL, ymax = asymp.UCL,
                     color = Salinity, group = Salinity)) +
  geom_point(position = position_dodge(width = 0.35)) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.35)) +
  geom_line(position = position_dodge(width = 0.35)) +
  facet_wrap(~ Soil, nrow = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Moisture level",
       y = "Predicted germination (%)",
       color = "Salinity") +
  theme_bw()
print(p_pred)
