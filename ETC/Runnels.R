
rm(list=ls())

library(dplyr)
library(tidyverse)

monitor1 <- read.csv("Reorganized_Runnels.csv")
monitor2 <- read.csv("InitialRunnels.csv")

monitor1 <- monitor1 %>%
  rename(Value_M1 = Value) 

monitor2 <- monitor2 %>%
  rename(Value_M2 = Value) 

monitor1 <- monitor1 %>%
  mutate(Value_M1 = as.numeric(Value_M1))

monitor2 <- monitor2 %>%
  mutate(Value_M2 = as.numeric(Value_M2))

names(monitor1)
names(monitor2)

comparison <- full_join(
  monitor1,
  monitor2,
  by = c(
    "Type.of.Plot",
    "Plot.Number",
    "Location.Type",
    "Location.Number",
    "Measurement.Category"
  )
) %>%
  mutate(Change = Value_M2 - Value_M1)

comparison %>%
  group_by(`Measurement.Category`) %>%
  summarise(
    Mean_M1 = mean(Value_M1, na.rm = TRUE),
    Mean_M2 = mean(Value_M2, na.rm = TRUE),
    Mean_Change = mean(Change, na.rm = TRUE),
    SD_Change = sd(Change, na.rm = TRUE),
    .groups = "drop"
  )
