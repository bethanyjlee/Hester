## using SeedSoaking.csv since that has daily change/new daily germ

##separating potting soil from hester mimic soil
seeds<-read.csv('SeedSoaking.csv')

##
 #install.packages("drcSeedGerm")
library(drcSeedGerm)
 #install.packages("GerminaR")
library(GerminaR)
library(tidyr)
library(dplyr)
library(knitr)

#overall germ analysis 
seed<-seeds %>%
  mutate(across(c(Soak.Time,Soak.Temp,Salinity.Soak,Soil.Type),as.factor))

gsm<-ger_summary(SeedN = "seeds", evalName = "D", data=seed)
gsm


#### creating Overall Data frames and graphs####
##data frame with percentage or relative germination in time by Soak.Temp
gtemp<-ger_intime(Factor = "Soak.Temp",
                  SeedN = "seeds",
                  evalName = "D",
                  method = "percentage",
                  data=seed)

tempgraph <- gtemp %>% 
  fplot(data = .
        , type = "line"
        , x = "evaluation"
        , y = "mean"
        , group = "Soak.Temp"
        , ylimits = c(0, 50, 10)
        , ylab = "Germination ('%')"
        , xlab = "Day"
        , glab = "Soaking Temperature (C)"
        , color = T
        , error = "ste"
  )

tempgraph
 ??fplot

##data frame with germination rate in time by soaking time
gtime<-ger_intime(Factor = "Soak.Time",
                  SeedN = "seeds",
                  evalName = "D",
                  method = "percentage",
                  data=seed)
gtime %>%
  head(10) %>%
  kable(caption= "Cumulative Germination by Soaking Time Factor")

timegraph <- gtime %>% 
  fplot(data = .
        , type = "line"
        , x = "evaluation"
        , y = "mean"
        , group = "Soak.Time"
        , ylimits = c(0, 50, 10)
        , ylab = "Germination ('%')"
        , xlab = "Day"
        , glab = "Soaking Time"
        , color = T
        , error = "ste"
  )

timegraph

##data frame with percentage or relative germination in time by Salinity Soak
gsalin<-ger_intime(Factor = "Salinity.Soak",
                  SeedN = "seeds",
                  evalName = "D",
                  method = "percentage",
                  data=seed)
gsalin %>%
  head(10) %>%
  kable(caption= "Cumulative Germination by Salinity Soak Factor")

salingraph <- gsalin %>% 
  fplot(data = .
        , type = "line"
        , x = "evaluation"
        , y = "mean"
        , group = "Salinity.Soak"
        , ylimits = c(0, 50, 10)
        , ylab = "Germination ('%')"
        , xlab = "Day"
        , glab = "Salinity Soak (ppt)"
        , color = T
        , error = "ste"
  )

salingraph


#Creating germination rate graphs but subsetting from soils

query1<-seeds$Soil.Type=='Potting Soil'
index1<-which(query1)
Potting<-seeds[index1,]
Mimic<-seeds[-index1,]

#Potting soil Data frames and graphs
# germ analysis 
pottingseeds<-Potting %>%
  mutate(across(c(Soak.Time,Soak.Temp,Salinity.Soak,Soil.Type),as.factor))

##data frame with percentage or relative germination in time by Soak.Temp
gtemppot<-ger_intime(Factor = "Soak.Temp",
                  SeedN = "seeds",
                  evalName = "D",
                  method = "percentage",
                  data=pottingseeds)
gtemppot %>%
  head(10) %>%
  kable(caption= "Potting Soil Germination by Soak Temperature Factor")

tempgraphpot <- gtemppot %>% 
  fplot(data = .
        , type = "line"
        , x = "evaluation"
        , y = "mean"
        , group = "Soak.Temp"
        , ylimits = c(0, 70, 10)
        , ylab = "Germination ('%')"
        , xlab = "Day"
        , glab = "Soaking Temperature (C)"
        , color = T
        , error = "ste"
  )

tempgraphpot


##data frame with germination rate in time by soaking time
gtimepot<-ger_intime(Factor = "Soak.Time",
                  SeedN = "seeds",
                  evalName = "D",
                  method = "percentage",
                  data=pottingseeds)
gtimepot %>%
  head(10) %>%
  kable(caption= "Potting Soil Germination by Soaking Time Factor")

timegraphpot <- gtimepot %>% 
  fplot(data = .
        , type = "line"
        , x = "evaluation"
        , y = "mean"
        , group = "Soak.Time"
        , ylimits = c(0, 70, 10)
        , ylab = "Germination ('%')"
        , xlab = "Day"
        , glab = "Soaking Time"
        , color = T
        , error = "ste"
  )

timegraphpot

##data frame with percentage or relative germination in time by Salinity Soak
gsalinpot<-ger_intime(Factor = "Salinity.Soak",
                   SeedN = "seeds",
                   evalName = "D",
                   method = "percentage",
                   data=pottingseeds)
gsalinpot %>%
  head(10) %>%
  kable(caption= "Potting Soil Germination by Salinity Soak Factor")

salingraphpot <- gsalinpot %>% 
  fplot(data = .
        , type = "line"
        , x = "evaluation"
        , y = "mean"
        , group = "Salinity.Soak"
        , ylimits = c(0, 70, 10)
        , ylab = "Germination ('%')"
        , xlab = "Day"
        , glab = "Salinity Soak (ppt)"
        , color = T
        , error = "ste"
  )

salingraphpot


#Hester mimic Data frames and graphs
# germ analysis 
mimicseeds<-Mimic %>%
  mutate(across(c(Soak.Time,Soak.Temp,Salinity.Soak,Soil.Type),as.factor))

##data frame with percentage or relative germination in time by Soak.Temp
gtempmim<-ger_intime(Factor = "Soak.Temp",
                     SeedN = "seeds",
                     evalName = "D",
                     method = "percentage",
                     data=mimicseeds)
gtempmim %>%
  head(10) %>%
  kable(caption= "Hester Mimic Soil Germination by Soak Temperature Factor")

tempgraphmim <- gtempmim %>% 
  fplot(data = .
        , type = "line"
        , x = "evaluation"
        , y = "mean"
        , group = "Soak.Temp"
        , ylimits = c(0, 40, 10)
        , ylab = "Germination ('%')"
        , xlab = "Day"
        , glab = "Soaking Temperature (C)"
        , color = T
        , error = "ste"
  )

tempgraphmim


##data frame with germination rate in time by soaking time
gtimemim<-ger_intime(Factor = "Soak.Time",
                     SeedN = "seeds",
                     evalName = "D",
                     method = "percentage",
                     data=mimicseeds)
gtimemim %>%
  head(10) %>%
  kable(caption= "Hester Mimic Soil Germination by Soaking Time Factor")

timegraphmim <- gtimemim %>% 
  fplot(data = .
        , type = "line"
        , x = "evaluation"
        , y = "mean"
        , group = "Soak.Time"
        , ylimits = c(0, 40, 10)
        , ylab = "Germination ('%')"
        , xlab = "Day"
        , glab = "Soaking Time"
        , color = T
        , error = "ste"
  )

timegraphmim

##data frame with percentage or relative germination in time by Salinity Soak
gsalinmim<-ger_intime(Factor = "Salinity.Soak",
                      SeedN = "seeds",
                      evalName = "D",
                      method = "percentage",
                      data=mimicseeds)
gsalinmim %>%
  head(10) %>%
  kable(caption= "Hester Mimic Soil Germination by Salinity Soak Factor")

salingraphmim <- gsalinmim %>% 
  fplot(data = .
        , type = "line"
        , x = "evaluation"
        , y = "mean"
        , group = "Salinity.Soak"
        , ylimits = c(0, 40, 10)
        , ylab = "Germination ('%')"
        , xlab = "Day"
        , glab = "Soaking Salinity"
        , color = T
        , error = "ste"
  )

salingraphmim


#### stats ####
#Check ANOVA assumptions
seedaov<-aov(mgr~Soak.Time+Soak.Temp+Salinity.Soak+Soil.Type,data=gsm)
summary(seedaov)
plot(seedaov)

#residuals vs fitted plot test for assumption 1 and 3 and were normal 
#homogeneity of variances

#Normal Q-Q graph shows left skewed tails
#normality of residuals 


#variable selections with no interactions 
DATA<-select(seed, contains("D"),-seeds,-ID.Number)

initial<-lm(mgr~1,data=gsm)
select<-stepAIC(initial,direction='forward',
                scope=list(upper= ~Soak.Time+Soak.Temp+Salinity.Soak+Soil.Type))
back<-stepAIC(seedaov,direction='backward',
              scope=list(upper= ~Soak.Time+Soak.Temp+Salinity.Soak+Soil.Type))
both<-stepAIC(seedaov,direction='both',
              scope=list(upper= ~Soak.Time+Soak.Temp+Salinity.Soak+Soil.Type))

#correct variables with no interactions 
seedaov2<-aov(mgr~Soak.Time+Salinity.Soak+Soil.Type,data=gsm)
summary(seedaov2)
plot(seedaov2)

#determine interactions 
intseed<-aov(mgr~Soak.Time*Salinity.Soak*Soil.Type,data=gsm)
summary(intseed)
plot(intseed)

#initially shows no interactions

#subset factors to make an interaction plot

germ<-gsm$mgr
time<-factor(gsm$Soak.Time)
saline<-factor(gsm$Salinity.Soak)
soil<-factor(gsm$Soil.Type)
temp<-factor(gsm$Soak.Temp)

#flat line in an interaction plot means no interaction 
#crossed lines + sig means interaction 
interaction.plot()

#
plot_grid(salingraphmim,salingraphpot, labels=c("Hester Mimic","Potting Soil"),)

p <- ggplot(seeds, aes(x=as.Date(, format = "Day"), y=Seedlings, color=Soil.Type, linetype = Salinity)) + 
  theme(
    panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(colour = "gray30"),
    axis.text.y = element_text(colour = "gray30"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size =20),
    axis.title.y = element_text(size =18)) +
  scale_colour_manual(values = c("tomato2", "yellow3","skyblue4","green4","darkslategray","grey60", "red", "orange"))+
  scale_y_continuous(limits = c(0, 25))+
  scale_shape_manual(values=c(0, 1, 2))+
  scale_linetype_manual(values=c("solid", "dotted", "dashed"))

#scale_colour_gradient(colours = c("yellow3", "turquoise3", "blue" ))

