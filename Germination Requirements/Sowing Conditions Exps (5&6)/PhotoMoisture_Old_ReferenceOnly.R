######  Zeanna Experiment - Moisture and Submersion
##### Analysis by Bethany Done on 8/19/25
####


####import data####
data<-read.csv('Zeanew.csv')

####load library####
library(ggplot2)
library(dplyr)
library(tidyr)
library(GerminaR)
library(brglm2)
library(emmeans)
library(car)

####pull columns needed and make long format####
#using max germ 
long<- data %>%
  pivot_longer(
    cols = starts_with("D"),
    names_to= "Date")

long<- long %>%
  mutate(Date= as.numeric(gsub("D","", Date)))

long<-long %>%
  mutate(value = as.double(value))

dataready<- long %>%
  mutate(across(c(Moisture_level_QL,Photo_period,exposure),as.factor))


#percentages
dataready1<-long %>%
  group_by(Moisture_level_QL,Photo_period,exposure) %>%
  summarise(germs = MaxGerm/Total)


##### determine interaction between moisture levels, photoperiod and exposure

interaction.plot(dataready1$Moisture_level_QL,dataready1$Photo_period,dataready1$germs)
###plot does not cross = no interaction 

interaction.plot(dataready1$Moisture_level_QL,dataready1$exposure,dataready1$germs)
#plot does not cross 

interaction.plot(dataready1$Photo_period,dataready1$exposure,dataready1$germs)
#plot does not cross



##lets see the differences in aic 

# Base (no interactions)
m0 <- glm(cbind(MaxGerm, Total - MaxGerm) ~ 
            Moisture_level_QL + Photo_period + exposure,
          family = binomial("logit"), data = dataready)

# Two-way interactions
m1 <- glm(cbind(MaxGerm, Total - MaxGerm) ~ 
            Moisture_level_QL * Photo_period + exposure,
          family = binomial("logit"), data = dataready)

m2 <- glm(cbind(MaxGerm, Total - MaxGerm) ~ 
            Moisture_level_QL * exposure + Photo_period,
          family = binomial("logit"), data = dataready)

m3 <- glm(cbind(MaxGerm, Total - MaxGerm) ~ 
            Photo_period * exposure + Moisture_level_QL,
          family = binomial("logit"), data = dataready)

# Full three-way interaction
m_full <- glm(cbind(MaxGerm, Total - MaxGerm) ~ 
                Moisture_level_QL * Photo_period * exposure,
              family = binomial("logit"), data = dataready)

# AIC comparison
AIC(m0, m1, m2, m3, m_full)

# Likelihood ratio tests
anova(m0, m1, test = "Chisq")
anova(m0, m2, test = "Chisq")
anova(m0, m3, test = "Chisq")
anova(m0, m_full, test = "Chisq")


#exact same AIC and same residual differene
deviance(m0) / df.residual(m0)
#### ahh this is super high, okay i guess we are running quasi again


#### stat analysis####
m0_quasi <- glm(cbind(MaxGerm, Total - MaxGerm) ~ 
                  Moisture_level_QL + Photo_period + exposure,
                family = quasibinomial("logit"), data = dataready)

summary(m0_quasi)

# overall tests

Anova(m0_quasi, type = 2)

# predicted probabilities

emm_m0_quasi <- emmeans(m0_quasi, ~ Moisture_level_QL + Photo_period + exposure,
                        type = "response")
emm_m0_quasi



####plot data ####

ggplot(data = dataready1, aes(x=Moisture_level_QL, y=germs, fill = Photo_period)) +
  geom_boxplot() +
  labs(y= "Germination Proportion")
