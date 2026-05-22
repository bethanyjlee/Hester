#### Chapter 2 Viability of Hester and Natural Sites
### elevation, succulence, site history

#### Old - prior to 9/26/25 - Analysis for CERF #####

setwd("~/R data/Chp2SeedSource/Exp 1 - Viability")
alldata<-read.csv('Viability.csv')

hester<-#need to filter all data to Site=HesterP1 and P2

hester %>%
  ggplot(aes(x=Elevation, y=Rate,)) +
  geom_point() +
  labs(title = "All Seed Viability by Elevation") +
 # facet_grid(~Site) +
  labs(x='Elevation', y='Percentage of Seeds Germinated') +
#  stat_summary(fun = "mean", color = 'red') + 
  theme_bw()

hestersum<-hester %>%
  group_by(Site,Elevation,Transect,Plot)%>%
  summarise(germs = mean(Rate))


# Fit a GLM (without random effects)
glm_model <- glm(cbind(Seeds, 50 - Seeds) ~ Elevation + Site, 
                 data = hester, 
                 family = binomial)

summary(glm_model)

glm_model2 <- glm(Seeds ~ Elevation*Site, data = hester, family = binomial)
summary(glm_model)

glm_model <- glm(cbind(Seeds, 50 - Seeds) ~ Elevation + Site + Succulence, 
                 data = hester, 
                 family = binomial)

# Fit a GLMM with plant random effects# Fit a GLMM with plant rsucculenceandom effects
library(lme4)
glmm_model <- glmer(Rate ~ Elevation + (1|IDTag), data = hester, family = binomial)

summary(glmm_model)

# Fit a GLMM with site random effects
glmm_model2 <- glm(Rate ~ Elevation + (1|Site), data = hester, family = binomial)
summary(glmm_model2)

# Compare models using ANOVA
anova(glm_model, glm_model2, test = "Chisq")


# Create a plot
ggplot(allseeds, aes(x = Elevation, y = Rate, color = SiteType)) +
  geom_point() +
  stat_smooth(method = "glm", method.args = list(family = "poisson"), se = TRUE) +
  theme_minimal() +
  labs(title = "Seed Quality of Salicornia by Elevation",
       x = "Elevation", y = "Germination Rate")


## Succulence Plots ##

succulence <- allseeds %>%
  filter(Site %in% c("Azevedo", "Hudson"))

               
glm_modelALL <- glm(cbind(Seeds, 50 - Seeds) ~ Elevation + Site + Succulence, data = allseeds,
                 family = binomial)

summary(glm_modelALL)


ggplot(succulence, aes(x = Elevation, y = Rate, color = Succulence, shape=Site, size = 2 )) +
  geom_point() +
  theme_minimal() +
  labs(title = "Germination Between Seed Succulence by Elevation",
       x = "Elevation", y = "Germination Rate")

allseeds$Sites <- factor(allseeds$Site, eallseeds$Sites <- factor(allseeds$Site, levels =c('Azevedo ','Hudson','OSRC','Yampah','HesterP1','HesterP2'))
allseeds$SucSite <- factor(allseeds$Site, levels =c('Azevedo ','Hudson'))

ggplot(allseeds, aes(x=Sites, y= Rate, color = Succulence)) +
  geom_boxplot() +
  labs(title = "Germination Between Seed Succulence by Site",
       x = "Site", y = "Germination Rate")


## Convert Elevation to Z scores ###

allseedsgerm$ZElev<-scale(allseedsgerm$Elevation, center = TRUE, scale = TRUE)
allseeds$ZElev<-scale(allseeds$Elevation, center = TRUE, scale = TRUE)

ggplot(allseedsgerm, aes(x = ZElev, y = germs, color = Succulence, shape=Site, size = 2 )) +
  geom_point() +
  theme_minimal() +
  labs(title = "Germination Between Seed Succulence by Elevation",
       x = "Elevation (meters)", y = "Germination Rate")

ggplot(allseeds, aes(x = ZElev, y = Rate, color = SiteType)) +
  geom_point()  +
  theme_minimal() +
  labs(title = "Seed Quality of Salicornia by Elevation",
       x = "Z-Score of Elevation ", y = "Germination Rate")

###### Background for Low, Medium and High Marsh ####

rects <- data.frame(
  xstart = c(-1.7, -0.5667, 0.5667),
  xend = c(-0.5667, 0.5667, 1.7),
  col = factor(c("Frequent", "King Tide Events Only", "Rarely to Never"), 
               levels = c("Rarely to Never", "King Tide Events Only","Frequent")))


ggplot(allseeds, aes(x = ZElev, y = Rate, color = SiteType)) +
  geom_point()  +
  theme_minimal() +
  labs(title = "Seed Viability of Salicornia pacifica by Elevation",
       x = "Z* for Elevation ", y = "Germination Rate") +
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col), 
            alpha = 0.2, inherit.aes = FALSE) +
  scale_fill_viridis_d(name = "Inundation Frequency", option = "C")


### parse out Natural vs Restored Marsh Graphs ###
  
ggplot(subset(allseeds, SiteType == "Natural"), aes(x = Elevation, y = Rate)) +
  geom_point(color = "red") +  
  labs(title = "Natural Sites - Seed Quality by Elevation",
       x = "Elevation (meters)", y = "Germination Rate") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "darkred", fill = "pink") +
  scale_y_continuous(limits = c(0,0.75)) +
  scale_x_continuous(limits = c(1.5,2.3)) +
  theme_minimal() 

ggplot(subset(allseeds, SiteType == "Restoration"), aes(x = Elevation, y = Rate)) +
  geom_point(color = "blue") +  
  labs(title = "Hester Marsh - Seed Quality by Elevation",
       x = "Elevation (meters)", y = "Germination Rate") +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), color = "darkred", fill = "pink") +
  scale_y_continuous(limits = c(0,0.75)) +
  scale_x_continuous(limits = c(1.5,2.3)) +
  theme_minimal() 


### phase one vs phase two ###

ggplot(subset(allseeds, SiteType == "Restoration"), aes(x = Elevation, y = Rate)) +
  geom_point(aes(color = Sites)) +  
  labs(title = "Hester Marsh - Seed Germination by Elevation",
       x = "Elevation (meters)", y = "Germination Rate") +
  geom_smooth(method = "glm", formula = y ~ poly(x,2), aes(color = Sites)) +
  scale_y_continuous(limits = c(0,0.75)) +
  scale_x_continuous(limits = c(1.5,2.2)) +
  theme_minimal() 


library(ggplot2)
library(dplyr)

# Define the background panels
rects <- data.frame(
  xstart = c(-1.7, -0.5667, 0.5667),
  xend = c(-0.5667, 0.5667, 1.7),
  col = factor(c("Frequent", "King Tide Events Only", "Rarely to Never"), 
               levels = c("Rarely to Never", "King Tide Events Only","Frequent"))
)

# Assign each ZElev to a panel category
allseeds <- allseeds %>%
  mutate(col = case_when(
    ZElev >= -1.7 & ZElev < -0.5667 ~ "Frequent",
    ZElev >= -0.5667 & ZElev < 0.5667 ~ "King Tide Events Only",
    ZElev >= 0.5667 & ZElev <= 1.7 ~ "Rarely to Never",
    TRUE ~ NA_character_
  )) %>%
  mutate(col = factor(col, levels = c("Rarely to Never", "King Tide Events Only","Frequent")))

# Calculate the mean germination rate per SiteType per panel
mean_rates <- allseeds %>%
  group_by(SiteType, col) %>%
  summarise(mean_rate = mean(Rate, na.rm = TRUE), .groups = "drop")

# Merge mean_rates with rects to get the midpoint of each category
mean_rates <- mean_rates %>%
  left_join(rects, by = "col") %>%
  mutate(mid_x = (xstart + xend) / 2)

# Plot with the averages added
ggplot(allseeds, aes(x = ZElev, y = Rate, color = SiteType)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Seed Viability of Salicornia pacifica by Elevation",
       x = "Z* for Elevation", y = "Germination Rate") +
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col), 
            alpha = 0.2, inherit.aes = FALSE) +
  scale_fill_viridis_d(name = "Inundation Frequency", option = "C") +
  geom_point(data = mean_rates, aes(x = mid_x, y = mean_rate, color = SiteType), 
             shape = 21, fill = "black", size = 6, stroke = 1, inherit.aes = FALSE)


#### New - 9/26/25 - Analysis for Manuscript ####


setwd("~/R data/Chp2SeedSource/Exp 1 - Viability")

#### Packages ####
library(dplyr)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(performance)
library(car)
library(emmeans)



# ---- Data ----
viab <- read.csv("Viability.csv") %>%
  mutate(
    Site       = factor(SiteID),
    SiteType   = factor(SiteType),       # Natural vs Restoration
 #   Succulence = factor(Succulence),     # Brown vs Green
    TidalCat   = factor(TidalCat),
    Fail       = Total - Seeds
  ) %>% filter(Total > 0)

data_summary <- viab %>%
# group_by(SiteID) %>%
 group_by(TidalCat) %>%
 # group_by(SiteType) %>%
  summarise(
    n        = n(),
    mean_germ = mean(Rate, na.rm = TRUE),
    sd_germ   = sd(Rate, na.rm = TRUE),
    se_germ   = sd_germ / sqrt(n)) %>%
  ungroup()

print(data_summary)


# 1) Base candidate models

##### Model Determination 1: SiteType - restored vs natural marshes #####

# Binomial GLM with Site fixed
m_type_null <- glm(cbind(Seeds, Fail) ~ 1,
                   data = viab, family = binomial)

m_type_main <- glm(cbind(Seeds, Fail) ~ SiteType + Elevation,
                   data = viab, family = binomial)

m_type_inter <- glm(cbind(Seeds, Fail) ~ SiteType * Elevation,
                    data = viab, family = binomial)

# Binomial GLMM with Site Randomized
m_type_glmm_null <- glmer(cbind(Seeds, Fail) ~ 1 + (1|Site),
                          data = viab, family = binomial,
                          control = glmerControl(optimizer = "bobyqa"))

m_type_glmm_main <- glmer(cbind(Seeds, Fail) ~ SiteType + Elevation + (1|Site),
                          data = viab, family = binomial,
                          control = glmerControl(optimizer = "bobyqa"))

m_type_glmm_inter <- glmer(cbind(Seeds, Fail) ~ SiteType * Elevation + (1|Site),
                           data = viab, family = binomial,
                           control = glmerControl(optimizer = "bobyqa"))




##### Model Determination 2: TidalCat - restored, tidal, or muted marshes #####

# Binomial GLM with Site fixed
m_tidal_null <- glm(cbind(Seeds, Fail) ~ 1,
                    data = viab, family = binomial)

m_tidal_main <- glm(cbind(Seeds, Fail) ~ TidalCat + Elevation,
                    data = viab, family = binomial)

m_tidal_inter <- glm(cbind(Seeds, Fail) ~ TidalCat * Elevation,
                     data = viab, family = binomial)


# Binomial GLMM 
m_tidal_glmm_null <- glmer(cbind(Seeds, Fail) ~ 1 + (1|Site),
                           data = viab, family = binomial,
                           control = glmerControl(optimizer = "bobyqa"))

m_tidal_glmm_main <- glmer(cbind(Seeds, Fail) ~ TidalCat + Elevation + (1|Site),
                           data = viab, family = binomial,
                           control = glmerControl(optimizer = "bobyqa"))

m_tidal_glmm_inter <- glmer(cbind(Seeds, Fail) ~ TidalCat * Elevation + (1|Site),
                            data = viab, family = binomial,
                            control = glmerControl(optimizer = "bobyqa"))


#### Model Determination 3: Site Only #####
m_site_main <- glm(cbind(Seeds, Fail) ~ Site + Elevation,
                   data = viab, family = binomial)

m_site_inter <- glm(cbind(Seeds, Fail) ~ Site * Elevation,
                    data = viab, family = binomial)


#### Model Selection ####
AIC(m_type_null, m_type_main, m_type_inter,
  m_type_glmm_null, m_type_glmm_main, m_type_glmm_inter,
  m_tidal_null, m_tidal_main, m_tidal_inter,
  m_tidal_glmm_null, m_tidal_glmm_main, m_tidal_glmm_inter)


# 2) Compare model fit

AIC(m_glm_site, m_glm_type, m_glmm_site, m_glm_full)

# Overdispersion checks

check_overdispersion(m_type_glmm_inter) ##around 2 - lets try beta binomial 

#### selected model ####
library(glmmTMB)

m_tidal_bb <- glmmTMB(
  cbind(Seeds, Fail) ~ TidalCat * Elevation + (1|Site),
  data = viab,
  family = betabinomial())

summary(m_tidal_bb)
performance::check_overdispersion(m_tidal_bb)

emmeans(m_tidal_bb, pairwise ~ TidalCat, type = "response")


##### 5) Posthoc summaries (EMMs) #####
emm_site <- emmeans(m_tidal_bb, ~ TidalCat * Elevation, type = "response")
print(emm_site)
print(pairs(emm_site, adjust = "tukey"))

#### making graphs ####

### Prediction grid ###
pred_grid <- expand.grid(
  Elevation = seq(min(viab$Elevation, na.rm = TRUE),
                  max(viab$Elevation, na.rm = TRUE),
                  length.out = 200),
  TidalCat = levels(viab$TidalCat))

pred_link <- predict(
  m_tidal_bb,
  newdata = pred_grid,
  type = "link",
  se.fit = TRUE,
  re.form = NA)

pred_grid <- pred_grid %>%
  mutate(
    fit_link = pred_link$fit,
    se_link  = pred_link$se.fit,
    fit      = plogis(fit_link),
    lwr      = plogis(fit_link - 1.96 * se_link),
    upr      = plogis(fit_link + 1.96 * se_link))


#### Plot 1 ####
### first make R2 and eqns ####
#### Extract equations + R2 from beta-binomial GLMM ####

library(performance)
#### Equations + marginal R2 ####

# Separate datasets
viab_tidal <- subset(viab, TidalCat == "Tidal")
viab_restored <- subset(viab, TidalCat == "Restored")

# Separate models
m_tidal <- glmmTMB(
  cbind(Seeds, Total - Seeds) ~ Elevation,
  family = betabinomial(link = "logit"),
  data = viab_tidal
)

m_restored <- glmmTMB(
  cbind(Seeds, Total - Seeds) ~ Elevation,
  family = betabinomial(link = "logit"),
  data = viab_restored
)

# Predicted values
pred_tidal <- predict(m_tidal, type = "response")
pred_restored <- predict(m_restored, type = "response")

# Pseudo R2
r2_tidal <- round(cor(viab_tidal$Seeds / viab_tidal$Total,
                      pred_tidal)^2, 3)

r2_restored <- round(cor(viab_restored$Seeds / viab_restored$Total,
                         pred_restored)^2, 3)

# Labels
eq_labels <- data.frame(
  TidalCat = c("Tidal", "Restored"),
  
  label = c(
    paste0(
      "Tidal: y = ",
      round(b0_tidal,3),
      ifelse(b1_tidal >= 0, " + ", " - "),
      abs(round(b1_tidal,3)),
      "x\nRm² = ", r2_tidal
    ),
    
    paste0(
      "Restored: y = ",
      round(b0_restored,3),
      ifelse(b1_restored >= 0, " + ", " - "),
      abs(round(b1_restored,3)),
      "x\nRm² = ", r2_restored
    )
  ),
  
  x = min(viab$Elevation) + 0.02,
  y = c(0.97, 0.88)
)


group_colors<- c(Tidal = "#1E88E5", Muted = '#D81B60', Restored = '#FFC107')

p_main <- ggplot() +
  geom_point(data = viab,
             aes(Elevation, Seeds/Total, color = TidalCat),
             alpha = 0.4, size = 2) +
 # geom_ribbon(data = pred_grid, aes(Elevation, ymin = lwr, ymax = upr, fill = TidalCat), alpha = 0.2) +
  geom_line(data = pred_grid,
            aes(Elevation, y = fit, color = TidalCat),
            linewidth = 1.2) +
  geom_text(data = eq_labels,aes(x = x,y = y,label = label,color = TidalCat),
            hjust = 0,size = 4,show.legend = FALSE)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0.02)) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 12)) +
  labs(x = "Elevation (m NAVD88)",
       y = "Proportion viable",
       color = "Tidal category",
       fill  = "Tidal category")

p_main
ggsave(filename = "Figure2.tif", dpi = 300, path = "Figures")

#### PLOT 1 but faceted ####

site_lookup <- viab %>%
  distinct(TidalCat, SiteType)

pred_grid <- pred_grid %>%
  left_join(site_lookup, by = "TidalCat")

neworder<-c("Tidal","Muted","Restored")
viab<-arrange(transform(viab, TidalCat = factor(TidalCat, levels = neworder)), TidalCat)

p_facet <- ggplot() +
  geom_point(data = viab,
             aes(Elevation, Seeds/Total, color = TidalCat),
             alpha = 0.4, size = 2) +
  geom_ribbon(data = pred_grid,
              aes(Elevation, ymin = lwr, ymax = upr, fill = TidalCat),
              alpha = 0.2) +
  geom_line(data = pred_grid,
            aes(Elevation, y = fit, color = TidalCat),
            linewidth = 1.2) +
  facet_wrap(~ TidalCat) +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_continuous(limits = c(0,0.8), expand = c(0,0.02)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 12)) +
  labs(title = "Seed viability across elevation by site type",
       x = "Elevation (m NAVD88)",
       y = "Proportion viable",
       color = "Tidal category")

p_facet


#### Supplemental PLOT 3 — Individual seed source / SiteID ####

# set site order
viab$Site <- factor(
  viab$Site,
  levels = c(
    "Restored",
    "Tidal1",
    "Tidal2",
    "Tidal3",
    "Tidal6"))

p_site <- ggplot(
  viab,
  aes(x = Elevation,
      y = Seeds / Total,
      color = Site)) +
  geom_point(alpha = 0.6, size = 2.5) +
  scale_y_continuous(
    limits = c(0, 1),expand = c(0, 0.02)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
    text = element_text(size = 12)) +
  labs( x = "Elevation (m NAVD88)",
    y = "Proportion viable")

p_site

ggsave(filename = "SupplementalFigure3.tif", dpi = 300, path = "Figures")


