
survseeds<-read.csv('Germination.csv')


library(GerminaR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)


str(survseeds)

  theme_bw()
long<- survseeds %>%
  pivot_longer(
    cols = starts_with("D"),
    names_to= "Date")

long<- long %>%
  mutate(Date= as.numeric(gsub("D","", Date)))

long<-long %>%
  mutate(value = as.numeric(value))

#counts
dataready1<- long %>%
  mutate(across(c(Soil.Type,Soak.Time,Soak.Temp,Salinity.Soak,Tray),as.factor))

#germination rate
plotsub<-dataready1 %>%
  group_by(Soil.Type,Soak.Time,Soak.Temp,Salinity.Soak,Date,Tray)%>%
  summarise(germ = value/seeds)

#average rates 
plotsub2<-plotsub %>%
  group_by(Date,Soil.Type,Soak.Time,Soak.Temp,Salinity.Soak)%>%
  summarise(germs = mean(germ))

plotsub2$SoakTime <- factor(plotsub2$Soak.Time, levels =c('3hr','1 day','7 days'))
plotsub2$SalSoak<- factor(plotsub2$Salinity.Soak, levels =c('Freshwater','15ppt','30ppt'))

plotsub2 %>%
  ggplot() +
  geom_line(aes(x=Date, y=germs, color = Soil.Type, linetype = Soak.Temp), linewidth = 1) +
  labs(title = "Seedling Count by Scarification") +
  labs(x='Days After Sown', y='Percentage of Seeds Alive') +
  facet_grid(SalSoak~SoakTime) +
  theme_classic()


#### ANOVA ####


soil<-plotsub$Soil.Type
time<-plotsub$Soak.Time
temp<-plotsub$Soak.Temp
sal<-plotsub$Salinity.Soak
germ<-plotsub$germ


model<-aov(germ~soil+temp+time*sal)
summary(model)
plot(model) ###this fits awful!!

## trying transformations ##

germsq<-sqrt(germ)

modelsq<-aov(germsq~soil+temp+time*sal)

plot(modelsq)

### make an interaction plot ###


interaction.plot(time, temp, germ)
#nope
interaction.plot(time,soil,germ)
#nope
interaction.plot(sal,temp,germ)
#nope
interaction.plot(sal,time,germ)
#SALSOAK AND SOAKTIME INTERACT
interaction.plot(sal,soil,germ)
#nope
interaction.plot(temp,soil,germ)
#nope
interaction.plot(soil,time,germ)
#nope
interaction.plot(tray,soil,germ)
#nope
interaction.plot(time,tray,germ)
#nope





##### NEW GERM ######

library(GerminaR)

new<-read.csv('seedSoaking.csv')

gsm<-ger_summary(SeedN = "seeds", evalName = "D", data=new)
summary(gsm)

#### by variable ####

query1<-new$Soil.Type=='Potting Soil'
index1<-which(query1)
Potting<-new[index1,]
Mimic<-new[-index1,]

hester<-ger_summary(SeedN = "seeds", evalName = "D", data=Mimic)
potting<-ger_summary(SeedN = "seeds", evalName = "D", data=Potting)

query2<-new$Salinity.Soak=='Freshwater'
index2<-which(query2)
Fresh<-new[index2,]

query3<-new$Salinity.Soak=='30ppt'
index3<-which(query3)
HighS<-new[index3,]

query4<-new$Salinity.Soak=='15ppt'
index4<-which(query4)
MedS<-new[index4,]

fresh<-ger_summary(SeedN = "seeds", evalName = "D", data=Fresh)
saltH<-ger_summary(SeedN = "seeds", evalName = "D", data=HighS)
saltM<-ger_summary(SeedN = "seeds", evalName = "D", data=MedS)

query5<-new$Soak.Temp=='12 C'
index5<-which(query5)
lowte<-new[index5,]

query6<-new$Soak.Temp=='20 C'
index6<-which(query6)
highte<-new[index6,]

LowT<-ger_summary(SeedN = "seeds", evalName = "D", data=lowte)
HighT<-ger_summary(SeedN = "seeds", evalName = "D", data=highte)

query7<-new$Soak.Time=='3hr'
index7<-which(query7)
threehr<-new[index7,]

query8<-new$Soak.Time=='1 day'
index8<-which(query8)
oneda<-new[index8,]

query9<-new$Soak.Time=='7 days'
index9<-which(query9)
sevend<-new[index9,]

threehours<-ger_summary(SeedN = "seeds", evalName = "D", data=threehr)
oneday<-ger_summary(SeedN = "seeds", evalName = "D", data=oneda)
sevendays<-ger_summary(SeedN = "seeds", evalName = "D", data=sevend)

summary(hester)
summary(potting)
summary(saltH)
summary(saltM)
summary(fresh)
summary(LowT)
summary(HighT)
summary(threehours)
summary(oneday)
summary(sevendays)

#### ANOVA for mgt ####

av <- aov(mgt ~ Salinity.Soak*Soak.Time+Soak.Temp+Soil.Type, data = gsm)
summary(av)

# mean comparison test

mc_mgt <- ger_testcomp(aov = av, comp = c("Salinity.Soak", "Soak.Time"))


mgt <- mc_mgt$table %>% 
  fplot(data = .
        , type = "bar" 
        , x = "Soak.Time"
        , y = "mgt"
        , group = "Salinity.Soak"
        , ylab = "Mean germination time (days)"
        , xlab = "Soaking Time"
        , glab = "Salinity Soak"
        , sig = "sig"
        , error = "ste"
        , color = T)

mc_mgt2 <- ger_testcomp(aov = av, comp = c("Soak.Temp", "Soak.Time"))


mc_mgt2$table %>% 
  fplot(data = .
        , type = "bar" 
        , x = "Soak.Time"
        , y = "mgt"
        , group = "Soak.Temp"
        , ylab = "Mean germination time (days)"
        , xlab = "Salinity"
        , glab = "Moisture Levels"
        , sig = "sig"
        , error = "ste"
        , color = T)

mc_mgt3 <- ger_testcomp(aov = av, comp = c("Soak.Temp", "Salinity.Soak"))


mc_mgt3$table %>% 
  fplot(data = .
        , type = "bar" 
        , x = "Salinity.Soak"
        , y = "mgt"
        , group = "Soak.Temp"
        , ylab = "Mean germination time (days)"
        , xlab = "Salinity"
        , glab = "Moisture Levels"
        , sig = "sig"
        , error = "ste"
        , color = T)


## potting soil only!! ##

query1<-new$Soil.Type=='Potting Soil'
index1<-which(query1)
Potting<-new[index1,]
Mimic<-new[-index1,]

PSgsm<-ger_summary(SeedN = "seeds", evalName = "D", data=Potting)

PSav <- aov(mgt ~ Salinity.Soak*Soak.Time+Soak.Temp, data = PSgsm)
summary(PSav)

# mean comparison test

PSmc_mgt <- ger_testcomp(aov = PSav, comp = c("Salinity.Soak", "Soak.Time"))

# data result




mc_mgt$table %>% 
  kable(caption = "Mean germination time comparison")


PSmgt <- PSmc_mgt$table %>% 
  fplot(data = .
        , type = "bar" 
        , x = "Soak.Time"
        , y = "mgt"
        , group = "Salinity.Soak"
        , ylab = "Mean germination time (days)"
        , xlab = "Soaking Time"
        , glab = "Salinity Soak"
        , sig = "sig"
        , color = F)

PSmgt

?fplot

### hester soil time

HMgsm<-ger_summary(SeedN = "seeds", evalName = "D", data=Mimic)

HMav <- aov(mgt ~ Salinity.Soak*Soak.Time+Soak.Temp, data = HMgsm)
summary(HMav)

# mean comparison test

HMmc_mgt <- ger_testcomp(aov = HMav, comp = c("Salinity.Soak", "Soak.Time"))

# data result


HMmgt <- HMmc_mgt$table %>% 
  fplot(data = .
        , type = "bar" 
        , x = "Soak.Time"
        , y = "mgt"
        , group = "Salinity.Soak"
        , ylab = "Mean germination time (days)"
        , xlab = "Soaking Time"
        , glab = "Salinity Soak"
        , sig = "sig"
        , color = T)

HMmgt


#### ANOVA for grs ####

## all soils

avg <- aov(grs ~ Salinity.Soak*Soak.Time+Soak.Temp+Soil.Type, data = gsm)
summary(avg)

# mean comparison test

mc_grs <- ger_testcomp(aov = avg, comp = c("Salinity.Soak", "Soak.Time"))

mc_grs$table$SalinityF = factor(mc_grs$table$Salinity.Soak, levels = c('Freshwater','15ppt','30ppt'))
mc_grs$table$SoakF = factor(mc_grs$table$Soak.Time, levels = c('3hr','1 day','7 days'))


grc <- mc_grs$table %>% 
  fplot(data = .
        , type = "bar" 
        , x = "SoakF"
        , y = "grs"
        , group = "SalinityF"
        , ylab = "Number of germinated seeds"
        , xlab = "Soaking Time"
        , glab = "All Soils        Salinity Soak"
        , sig = "sig"
        , color = F)
grc

?fplot

## potting soil only!! ##

PSavg <- aov(grs ~ Salinity.Soak*Soak.Time+Soak.Temp, data = PSgsm)
summary(PSavg)
plot(PSavg)

# mean comparison test

PS_grs <- ger_testcomp(aov = PSavg, comp = c("Salinity.Soak", "Soak.Time"))


PS_grs$table$SalinityF = factor(PS_grs$table$Salinity.Soak, levels = c('Freshwater','15ppt','30ppt','Control'))
PS_grs$table$SoakF = factor(PS_grs$table$Soak.Time, levels = c('3hr','1 day','7 days','Control'))

table<-PS_grs$table

PSgrc <- PS_grs$table %>% 
  fplot(data = .
        , type = "bar" 
        , x = "SoakF"
        , y = "grs"
        , group = "SalinityF"
        , ylab = "Number of germinated seeds"
        , xlab = "Soaking Time"
        , glab = "Potting Soil            Salinity Soak"
        , sig = "sig"
        , color = F)

PSgrc

ggplot(data = table, aes(x=SoakF, y=grs, fill = SalinityF)) +
  geom_col(position = position_dodge2(preserve = 'single')) + 
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Seeds Germinated in High Nutrient, Potting Soil", x = "Soaking Time",
       y = "Number of Germinated Seeds", fill = "Salinity Soak")

#### relevant PS plot for MS ####

PSgsm$SoakF = factor(PSgsm$Salinity.Soak, levels = c('Freshwater','15ppt','30ppt','Control'))
PSgsm$TimeF = factor(PSgsm$Soak.Time, levels = c('3hr','1 day','7 days','Control'))

ggplot(data = PSgsm, aes(x = TimeF, y = grs, fill = SoakF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) +  # Bar width adjustment
  scale_fill_viridis_d(option = "C") +
  labs(title = "Seeds Germinated in High Nutrient, Potting Soil", 
       x = "Soaking Time", 
       y = "Number of Germinated Seeds", 
       fill = "Salinity Soak") +
  scale_y_continuous(limits = c(0,35))+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw()


PSMGT<-ggplot(data = PSgsm, aes(x = TimeF, y = mgt, fill = SoakF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) +  # Bar width adjustment
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Days to Germinate in High Nutrient, Potting Soil", 
       x = "Soaking Time", 
       y = "Days to Germinate", 
       fill = "Salinity Soak") +
  scale_y_continuous(limits = c(0,8))+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw()

### hester soil time

HMavg <- aov(grs ~ Salinity.Soak*Soak.Time+Soak.Temp, data = HMgsm)
summary(HMavg)

plot(HMavg)

# mean comparison test

HMmc_grs <- ger_testcomp(aov = HMavg, comp = c("Salinity.Soak", "Soak.Time"))

HMmc_grs$table$SalinityF = factor(HMmc_grs$table$Salinity.Soak, levels = c('Freshwater','15ppt','30ppt','Control'))
HMmc_grs$table$SoakF = factor(HMmc_grs$table$Soak.Time, levels = c('3hr','1 day','7 days','Control'))

table2<-HMmc_grs$table

HMgrs <- HMmc_grs$table %>% 
  fplot(data = .
        , type = "bar" 
        , x = "SoakF"
        , y = "grs"
        , group = "SalinityF"
        , ylab = "Number of germinated seeds"
        , xlab = "Soaking Time"
        , glab = "Hester Mimic           Salinity Soak"
        , sig = "sig"
        , color = F)

HMgrs

#### relevant HM plot for MS ####

ggplot(data = table2, aes(x=SoakF, y=grs, fill = SalinityF)) +
  geom_col(position = position_dodge2(preserve = 'single')) + 
  scale_fill_viridis_d(option = "C") +
  scale_y_continuous(limits = c(0,35))+
  labs(title = "Seeds Germinated in Low Nutrient, Clay Soil", x = "Soaking Time",
       y = "Number of Germinated Seeds", fill = "Salinity Soak")+
  theme_bw()


HMgsm$SoakF = factor(HMgsm$Salinity.Soak, levels = c('Freshwater','15ppt','30ppt','Control'))
HMgsm$TimeF = factor(HMgsm$Soak.Time, levels = c('3hr','1 day','7 days','Control'))

ggplot(data = HMgsm, aes(x = TimeF, y = grs, fill = SoakF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) +  # Bar width adjustment
  scale_fill_viridis_d(option = "C") +
  labs(title = "Seeds Germinated in Nutrient Poor, Hester Soil", 
       x = "Soaking Time", 
       y = "Number of Germinated Seeds", 
       fill = "Salinity Soak") +
  scale_y_continuous(limits = c(0,35))+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw()



HMMGT<-ggplot(data = HMgsm, aes(x = TimeF, y = mgt, fill = SoakF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) +  # Bar width adjustment
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Days to Germinate in Nutrient Poor, Hester Soil", 
       x = "Soaking Time", 
       y = "Days to Germinate", 
       fill = "Salinity Soak") +
  scale_y_continuous(limits = c(0,8))+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw()

library(ggplot2)


####### redone data - start here??? ######



library(ggplot2)
library(dplyr)
library(gridExtra)

# Step 1: Replace zeros with NA in germination columns
new_clean <- new %>%
  mutate(across(starts_with("D"), ~ ifelse(. == 0, NA, .)))

# Step 2: Recalculate summaries for Potting Soil and Hester Mimic soil
PSgsm_clean <- ger_summary(SeedN = "seeds", evalName = "D", data = new_clean %>% filter(Soil.Type == "Potting Soil"))
HMgsm_clean <- ger_summary(SeedN = "seeds", evalName = "D", data = new_clean %>% filter(Soil.Type == "Hester Mimic"))

# Step 3: Handle missing data (NA) and ensure proper column types
PSgsm_clean <- PSgsm_clean %>%
  filter(!is.na(mgt) & !is.na(grs)) %>%
  mutate(TimeF = as.factor(TimeF), SoakF = as.factor(SoakF))

HMgsm_clean <- HMgsm_clean %>%
  filter(!is.na(mgt) & !is.na(grs)) %>%
  mutate(TimeF = as.factor(TimeF), SoakF = as.factor(SoakF))

# Step 4: Set the black and white color palette
bw_palette <- c("black", "grey20", "grey60", "white")

# Plot 1: Mean Germination Time for Potting Soil
PSMGT <- ggplot(data = PSgsm_clean, aes(x = TimeF, y = mgt, fill = SoakF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) + 
  scale_fill_manual(values = bw_palette) +
  labs(title = "Days to Germinate in High Nutrient, Potting Soil", 
       x = "Soaking Time", 
       y = "Days to Germinate", 
       fill = "Salinity Soak") +
  scale_y_continuous(limits = c(0, 8)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw() + theme(legend.position = "none")

# Plot 2: Mean Germination Time for Hester Mimic Soil
HMMGT <- ggplot(data = HMgsm_clean, aes(x = TimeF, y = mgt, fill = SoakF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) + 
  scale_fill_manual(values = bw_palette) +
  labs(title = "Days to Germinate in Nutrient Poor, Hester Soil", 
       x = "Soaking Time", 
       y = "Days to Germinate", 
       fill = "Salinity Soak") +
  scale_y_continuous(limits = c(0, 8)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw() + theme(legend.position = "none")

# Plot 3: Mean Germinated Seeds for Potting Soil
PSgrc <- ggplot(data = PSgsm_clean, aes(x = TimeF, y = grs, fill = SoakF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) + 
  scale_fill_manual(values = bw_palette) +
  labs(title = "Number of Germinated Seeds in High Nutrient, Potting Soil", 
       x = "Soaking Time", 
       y = "Number of Germinated Seeds", 
       fill = "Salinity Soak") +
  scale_y_continuous(limits = c(0, 35)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw() + theme(legend.position = "none")

# Plot 4: Mean Germinated Seeds for Hester Mimic Soil
HMgrs <- ggplot(data = HMgsm_clean, aes(x = TimeF, y = grs, fill = SoakF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) + 
  scale_fill_manual(values = bw_palette) +
  labs(title = "Number of Germinated Seeds in Nutrient Poor, Hester Soil", 
       x = "Soaking Time", 
       y = "Number of Germinated Seeds", 
       fill = "Salinity Soak") +
  scale_y_continuous(limits = c(0, 35)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw() + theme(legend.position = "none")

# Step 5: Combine all the plots in one figure with one shared legend on the right
grid.arrange(PSMGT, HMMGT, PSgrc, HMgrs, ncol = 2,
             top = "Comparison of Germination Time and Seed Count", 
             bottom = "Soaking Time",
             left = "Mean Germination Time and Seed Count") 


##### summary data - pulling what we need from this ####

#  1) Summary data from GerminaR output

# gsm is the ger_summary() output
gsm <- gsm %>%
  mutate(
    Soil.Type     = factor(Soil.Type),
    Soak.Time     = factor(Soak.Time, levels = c("3hr","1 day","7 days","Control")),
    Salinity.Soak = factor(Salinity.Soak, levels = c("Freshwater","15ppt","30ppt","Control")),
    Soak.Temp     = factor(Soak.Temp)
  )

# Quick grouped summary (means ± SE) for reporting
sum_by_trt <- gsm %>%
  group_by(Soil.Type, Salinity.Soak, Soak.Time, Soak.Temp) %>%
  summarise(
    n         = n(),
    seedsN    = mean(seeds, na.rm = TRUE),
    grs_mean  = mean(grs, na.rm = TRUE),
    grs_se    = sd(grs, na.rm = TRUE) / sqrt(n),
    mgt_mean  = mean(mgt, na.rm = TRUE),
    mgt_se    = sd(mgt, na.rm = TRUE) / sqrt(n),
    .groups   = "drop"
  )

print(head(sum_by_trt), row.names = FALSE)

## ----------------------------
## 2) Model A: Germination proportion (binomial)
## ----------------------------
# Response: germinated counts vs. non-germinated
mod_GP <- glm(
  cbind(grs, seeds - grs) ~ Soil.Type + Soak.Temp + Soak.Time * Salinity.Soak,
  data = gsm,
  family = binomial
)

cat("\n=== Binomial model for germination proportion ===\n")
summary(mod_GP)

## ----------------------------
## 3) Model B: Mean Germination Time (Gamma)
## ----------------------------
# Response: mgt from GerminaR (positive + right-skewed)
mod_mgt_gamma <- glm(
  mgt ~ Soil.Type + Soak.Temp + Soak.Time * Salinity.Soak,
  data = gsm,
  family = Gamma(link = "log")
)

cat("\n=== Gamma(log) model for mean germination time ===\n")
summary(mod_mgt_gamma)
