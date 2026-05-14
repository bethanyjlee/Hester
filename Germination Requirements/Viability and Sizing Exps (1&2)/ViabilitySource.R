#### Experiment for seed source testing 

### Viability tests were done by the Baskin2 methodology and were done 
## every few months post October/November Harvest for a year

setwd("~/R data/Chp1PWSeeds/ViabilityANDSizing")

#### library ####

library(tidyverse) 
library(dplyr)
library(car)
library(ggpubr)
library(ggrepel)
library(grid)
library(multcompView)

#### Viability Tests ####



####June#####
viability<-read.csv('JuneViability.csv')



viability<-viability[!grepl('Late ',viability$Harvest),]

viability %>%
  ggplot(aes(x=Seed.Source, y= Germ))+
  geom_boxplot()+
  labs(title = "Viability Results 8 Months After Harvest", x= "Seed Source", y= "Germination Percentage")+
  scale_y_continuous(labels = scales::percent, limits = c(0,0.6)) +
  theme_bw()


avg<- viability %>%
  group_by(Seed.Source,Harvest) %>%
  summarize(Germs = mean(Germ)) 
 
#### old Tukeys HSD test ####

oneway<-aov(Germ ~ Seed.Source, data = viability)

summary(oneway)

tukey<-TukeyHSD(oneway,"Seed.Source")

cld<-multcompLetters4(oneway,tukey)

cld<-as.data.frame.list(cld$Seed.Source)

dt<- group_by(viability, Seed.Source) %>%
  summarise(size = mean(viability$Germ), sd = sd(viability$Germ)) 

dt$cld<-cld$Letters
print(dt$cld)


?geom_text

#### new tukey hsd test ####


# Boxplot of initial viability results
viability %>%
  ggplot(aes(x = Seed.Source, y = Germ)) +
  geom_boxplot() +
  labs(title = "Initial Viability Results After Harvest", 
       x = "Seed Source", y = "Germination Percentage") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.6)) +
  theme_bw()

# Calculate means for annotation
dt <- viability %>%
  group_by(Seed.Source) %>%
  summarise(mean_germ = mean(Germ, na.rm = TRUE), 
            max_germ = max(Germ, na.rm = TRUE),  # Get max for annotation position
            .groups = "drop")



# One-way ANOVA
oneway <- aov(Germ ~ Seed.Source, data = viability)
summary(oneway)

sd(viability$Germ)

# Tukey HSD Test
tukey <- TukeyHSD(oneway, "Seed.Source")

# Extract significance letters
cld <- multcompLetters4(oneway, tukey)

# Convert letters into a dataframe
cld_df <- data.frame(Seed.Source = names(cld$Seed.Source$Letters), 
                     Letters = cld$Seed.Source$Letters)

# Merge letters with dt
dt <- left_join(dt, cld_df, by = "Seed.Source")

# Print letters to check
print(dt$Letters)

# **Plot with Tukey Letters**
ggplot(viability, aes(x = Seed.Source, y = Germ)) +
  geom_boxplot() +
  geom_text(data = dt, aes(x = Seed.Source, y = max_germ + 0.01, label = Letters), 
            size = 6, vjust = 0) +
  labs(title = "Initial Viability Results After Harvest", 
       x = "Seed Source", y = "Germination Percentage") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.65)) +
  theme_bw()



######November Harvest########


# Load data
viability2 <- read.csv('OctoberViability.csv')

# Remove 'Late' harvest data
viability2 <- viability2[!grepl('Late', viability2$Harvest),]


# Boxplot of viability one year after harvest
ggplot(viability2, aes(x = Seed.Source, y = Germ)) +
  geom_boxplot() +
  labs(title = "Viability One Year After Harvest", 
       x = "Seed Source", y = "Germination Percentage") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.6)) +
  theme_bw()

# Summarize data for Tukey annotation
dt2 <- viability2 %>%
  group_by(Seed.Source) %>%
  summarise(mean_germ = mean(Germ, na.rm = TRUE), 
            max_germ = max(Germ, na.rm = TRUE),  # Get max for annotation position
            .groups = "drop")

# One-way ANOVA
oneway2 <- aov(Germ ~ Seed.Source, data = viability2)
summary(oneway2)  # FIXED: Using `oneway2` instead of `oneway`

sd(viability2$Germ)

# Tukey HSD Test
tukey2 <- TukeyHSD(oneway2, "Seed.Source")

# Extract significance letters
cld2 <- multcompLetters4(oneway2, tukey2)

# Convert Tukey letters into a dataframe
cld_df2 <- data.frame(Seed.Source = names(cld2$Seed.Source$Letters),  # FIXED: Use correct object
            Letters = cld2$Seed.Source$Letters)  # FIXED: Reference `cld2`, not `cld`

# Merge Tukey letters with summary data
dt2 <- left_join(dt2, cld_df2, by = "Seed.Source")

# Print letters to check
print(dt2$Letters)

# **Final Boxplot with Tukey Letters**
ggplot(viability2, aes(x = Seed.Source, y = Germ)) +
  geom_boxplot() +
  geom_text(data = dt2, aes(x = Seed.Source, y = max_germ + 0.01, label = Letters),  # FIXED: Use `dt2`
            size = 6, vjust = 0) +
  labs(title = "Pickleweed Viability One Year After Harvest", 
       x = "Seed Source", y = "Germination Percentage") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.65)) +
  theme_bw()

####Harvest#####

early<-viability2[!grepl('Late',viability2$Harvest),]
late<-viability2[!grepl('Early',viability2$Harvest),]

filtered_data <- viability2 %>%
  filter(Seed.Source %in% c("Hester2", "MoroCojo"))


ggplot(filtered_data, aes(x = Seed.Source, y = Germ)) +
  geom_boxplot() +
  labs(title = "Viability One Year After Harvest", 
       x = "Seed Source", y = "Germination Percentage") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.65)) +
  facet_wrap(~Harvest) +
  theme_bw()


