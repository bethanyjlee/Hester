##### Rate to Size ####
sized<-read.csv('GermSized.csv')



mean(sized$Sized)

median(sized$Sized)

early<-sized[!grepl('Late',sized$Harvest),]
late<-sized[!grepl('Early',sized$Harvest),]
aggregate(Sized ~ Seed.Source, data = late, FUN = mean, na.rm = TRUE)

mean(late$Sized)

median(late$Sized)

early_subset <- sized[sized$Harvest == "Early" & sized$Seed.Source %in% c("Hester2", "MoroCojo"), ]
aggregate(Sized ~ Seed.Source, data = early_subset, FUN = mean, na.rm = TRUE)

late_subset <- sized[sized$Harvest == "Late" & sized$Seed.Source %in% c("Hester2", "MoroCojo"), ]
aggregate(Sized ~ Seed.Source, data = late_subset, FUN = mean, na.rm = TRUE)

#### subset both sources and run a t-test to compare mean sizes ####

moro_data <- sized[sized$Seed.Source == "MoroCojo", ]

shapiro.test(sized$Sized[sized$Seed.Source == "MoroCojo" & sized$Harvest == "Early"]) #not normal
wilcox.test(Sized ~ Harvest, data = moro_data)


## run a multiple regression ##
## check assumptions ##
lmGerm<-lm(Germ ~ Sized, data = early)
summary(lmGerm)

plot(lmGerm)


lmGerm2<-lm(Germ ~ AvgSize + Seed.Source, data = early)
summary(lmGerm2)
plot(lmGerm2)

vif(lmGerm2) ### small vif = run together

aovsize <- aov(Germ ~ AvgSize + Seed.Source, data = early)
summary(aovsize)

## run an interaction ##

lmGerm3<- aov(Germ ~ Seed.Source*Sized, data = early)
summary(lmGerm3)
read

#plot data#

rmse <- round(sqrt(mean(resid(lmGerm2)^2)), 2)
coefs <- unname(coef(lmGerm2))  # Remove names
b0 <- round(coefs[1], 2)
b1 <- round(coefs[2],2)
r2 <- round(summary(lmGerm2)$r.squared, 2)  # Round R² for better display

r2

#eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
##             r^2 == .(r2) * "," ~~ RMSE == .(rmse))

# eqn <- as.expression(bquote(italic(y) == .(b0) + .(b1) * italic(x) * "," ~~ 
##                           R^2 == .(r2) * "," ~~ RMSE == .(rmse)))

# eqn <- deparse(bquote(italic(y) == .(b0) + .(b1) * italic(x) * "," ~~ 
##                     R^2 == .(r2) * "," ~~ RMSE == .(rmse)))


eqn <- paste0("italic(y) == ", b0, 
              ifelse(b1 >= 0, " + ", " - "), abs(b1), " * italic(x) * ',' ~ ~",
              "R^2 == ", r2, " * ',' ~ ~",
              "RMSE == ", rmse)

eqn

early1<- early %>%
  group_by(Seed.Source,AvgSize,Replicate) %>%
  summarize(Germs = mean(Germ))


early1 %>%
  ggplot(aes(x= AvgSize, y=Germs)) +
  geom_point(aes(color = Seed.Source))+
  labs(title = "Relationship of Germination Rate to Size",
       x = "Average Size (mm)", y= "Germination Percentage")+
  geom_smooth(method = lm, se =FALSE)+
  #  annotate("text", x = 1.5, y = 0.5, label = "R² = 0.68")+
  annotate("text", x = 1.5, y = 0.55, label = eqn, parse = TRUE) +
  scale_y_continuous(labels = scales::percent, limits = c(0,0.6)) +
  theme(legend.position = "none") +
  theme_classic()


?annotate



# Boxplot of AvgSize across Seed.Source

ggplot(early, aes(x = Seed.Source, y = AvgSize)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Distribution of Seed Size by Seed Source")



#### GLM ####

### as the fata is proportion we will use binomial with weights

TotalSeeds <- 50

glm_model <- glm(cbind(Germ * TotalSeeds, (1 - Germ) * TotalSeeds) ~ AvgSize + Seed.Source, 
                 data = early, family = binomial)
summary(glm_model)



aggregate(Germ ~ Seed.Source + Harvest, data = sized, FUN = mean, na.rm = TRUE)


##### 8/11/25 analysis ######

library(tidyverse)

# Read data
sized <- read.csv("GermSized.csv")

#--- 1) Size differences among sites ---
size_aov <- aov(Sized ~ Seed.Source, data = sized)
summary(size_aov)
TukeyHSD(size_aov) # Post-hoc differences

#--- 2) Relationship between size and germination ---
# We'll assume Germ is proportion germinated (0–1) and each trial used TotalSeeds seeds
TotalSeeds <- 50  

# Model with site
glm_site <- glm(cbind(Germ * TotalSeeds, (1 - Germ) * TotalSeeds) ~ Sized + Seed.Source,
                data = sized, family = binomial)

# Model without site
glm_nosite <- glm(cbind(Germ * TotalSeeds, (1 - Germ) * TotalSeeds) ~ Sized,
                  data = sized, family = binomial)

# Compare models
AIC(glm_site, glm_nosite)
anova(glm_nosite, glm_site, test = "Chisq") # LRT to test if site improves fit

# Summaries
summary(glm_site)
summary(glm_nosite)

#--- 3) Scatterplot of size vs germination % ---
ggplot(sized, aes(x = Sized, y = Germ, color = Seed.Source)) +
  geom_point(size = 2) +
geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
  labs(x = "Seed Size (mm)",
       y = "Germination Percentage",
       color = "Site") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic()



