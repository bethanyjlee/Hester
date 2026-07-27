rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)

dat<-read.csv('Walled_UAV.csv')

boxplot(PercentVeg~Walled, data = dat)

hist(dat$PercentVeg[dat$Walled == "Yes"],
     main = "Walled",
     xlab = "% Vegetation")

hist(dat$PercentVeg[dat$Walled == "No"],
     main = "Walled",
     xlab = "% Vegetation")

shapiro.test(dat$PercentVeg[dat$Walled == "No"])

shapiro.test(dat$PercentVeg[dat$Walled == "Yes"])

var.test(PercentVeg ~ Walled, data = dat) ### p value of 0.3652 = variances are simialar

qqnorm(dat$PercentVeg[dat$Walled == "No"])
qqline(dat$PercentVeg[dat$Walled == "No"])

qqnorm(dat$PercentVeg[dat$Walled == "Yes"])
qqline(dat$PercentVeg[dat$Walled == "Yes"])

#### normality is clearly violated, lets run wilxoc test

wilcox.test(PercentVeg ~ Walled, data = dat)

ggplot(dat, aes(x = Walled, y = PercentVeg, fill = Walled)) +
  geom_boxplot() + 
  geom_jitter() +
  labs(x = "Walled or Not", y = "Vegetation cover %") +
  theme_bw()+
  theme(legend.position = "none") +
  facet_wrap(~GARDENNUM)
