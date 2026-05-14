library(GerminaR)
library(ggplot2)
library(tidyr)
library(lubridate)
library(MASS)
library(car)
library(ggpubr)
library(knitr)

 
## before ####
pilot<-read.csv('rawdata.csv')

pilot1<-pilot %>%
  mutate(across(c(Soil,Compaction,Salinity,Moisture)))

gsm<-ger_summary(SeedN = "Seeds", evalName = "D", data=fulldata)
gsm


gsalinity<-ger_intime(Factor = "Salinity",
                  SeedN = "Seeds",
                  evalName = "D",
                  method = "percentage",
                  data=pilot1)


gsalgraph<-gsalinity %>% 
  fplot(data = .
        , type = "line"
        , x = "evaluation"
        , y = "mean"
        , group = "Salinity"
        , ylimits = c(0, 100, 10)
        , ylab = "Germination ('%')"
        , xlab = "Day"
        , glab = "Soaking Salinity (ppt)"
        , color = T
        , error = "ste"
  )

gsalgraph



query1<-pilot1$Soil=='PottingSoil'
index1<-which(query1)
Potting<-pilot1[index1,]
Mimic<-pilot1[-index1,]

## Potting soil  ####
gsalinitypot<-ger_intime(Factor = "Salinity",
                      SeedN = "Seeds",
                      evalName = "D",
                      method = "percentage",
                      data=Potting)


gsalgraphpot<-gsalinitypot %>% 
  fplot(data = .
        , type = "line"
        , x = "evaluation"
        , y = "mean"
        , group = "Salinity"
        , ylimits = c(0, 100, 10)
        , ylab = "Germination ('%')"
        , xlab = "Day"
        , glab = "Soaking Salinity (ppt)"
        , color = T
        , error = "ste"
  )

gsalgraphpot


##Hester Mimic Soil  ####
## Potting soil 
gsalinitymim<-ger_intime(Factor = "Salinity",
                         SeedN = "Seeds",
                         evalName = "D",
                         method = "percentage",
                         data=Mimic)


gsalgraphmim<-gsalinitymim %>% 
  fplot(data = .
        , type = "line"
        , x = "evaluation"
        , y = "mean"
        , group = "Salinity"
        , ylimits = c(0, 100, 10)
        , ylab = "Germination ('%')"
        , xlab = "Day"
        , glab = "Soaking Salinity (ppt)"
        , color = T
        , error = "ste"
  )

gsalgraphmim



## check anova ####

seedaov<-aov(mgr~Salinity+Moisture+Compaction+Soil,data=gsm)
summary(seedaov)
plot(seedaov)
trial1<-plot(seedaov,alpha=0.05,LSD.test("Soak.Temp",seeds))



## code from viz final ####

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

allav<- aov(mgt~Salinity*Moisture*Soil, data = gsm)
summary(allav)
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


bw_palette <- c("black", "grey20", "grey60", "grey80")


PSGRS<-ggplot(data = potting, aes(x = Salinity, y = grs, fill = MoistureF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) +  # Bar width adjustment
 # scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values = bw_palette)+
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

HMGRS <- ggplot(data = hester, aes(x = Salinity, y = grs, fill = MoistureF)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), width = 0.7) +
 # scale_fill_brewer(palette = "Set2")+
  scale_fill_manual(values = bw_palette)+
  labs(title = "Seeds Germinated in Nutrient Poor, Hester Soil", 
       x = "Moisture Addition to Soils", 
       y = "Number of Germinated Seeds", 
       fill = "Moisture Levels") +
  scale_y_continuous(limits = c(0,35)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8)) + 
  theme_bw()



HMGRS
PSGRS
