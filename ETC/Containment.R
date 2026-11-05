library(readxl)
library(dplyr)
library(ggplot2)
library(emmeans)

dat <- read.csv("Containment.csv") %>%
  mutate(
    Site = factor(Site),
    Elevation = factor(Elevation,
                       levels = c("LowLow", "Low", "MidLow", "Middle",
                                  "MidHigh", "High", "UpperHigh")),
    Containment = relevel(factor(Containment), ref = "Burlap"),
    AllSeedsEnd = as.integer(AllSeedsEnd),
    SeedsStart = as.integer(SeedsStart)
  )

# Order elevations from lowest to highest
elevation_order <- c(
  "LowLow", "Low", "MidLow", "Middle",
  "MidHigh", "High", "UpperHigh"
)

# Clean data and calculate containment percentage
plot_data <- dat %>%
  mutate(
    Elevation = factor(Elevation, levels = elevation_order),
    Containment = factor(
      Containment,
      levels = c("PeatPot", "Burlap", "WoodPulpBag", "CornFiberBag")
    ),
    Site = factor(Site),
    Containment. = if_else(
      SeedsStart > 0,
      AllSeedsEnd / SeedsStart * 100,
      NA_real_
    )
  )
 
# Consistent graph theme
my_theme <- theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 30, hjust = 1, color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  )

# GRAPH 1: Containment percentage by method
ggplot(
  plot_data,
  aes(x = Containment, y = Containment., fill = Containment)
) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.7, size = 2) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = "Containment Performance by Method",
    x = "Containment Method",
    y = "Seeds Retained (%)"
  ) +
  my_theme +
  theme(legend.position = "none")


# GRAPH 2: Containment percentage across elevations
ggplot(
  plot_data,
  aes(
    x = Elevation,
    y = Containment.,
    color = Containment,
    group = Containment
  )
) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.15
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
  labs(
    title = "Containment Performance Across Elevations",
    x = "Elevation",
    y = "Mean Seeds Retained (%)"
  ) +
  my_theme


# GRAPH 3: Containment percentage by site and method
ggplot(
  plot_data,
  aes(
    x = Site,
    y = Containment.,
    fill = Containment
  )
) +
  stat_summary(
    fun = mean,
    geom = "col",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = "Containment Performance by Site",
    x = "Site",
    y = "Mean Seeds Retained (%)"
  ) +
  my_theme


# GRAPH 4: Mean number of seeds remaining by elevation and method
ggplot(
  plot_data,
  aes(
    x = Elevation,
    y = AllSeedsEnd,
    color = Containment,
    group = Containment
  )
) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Mean Seeds Remaining Across Elevations",
    x = "Elevation",
    y = "Mean Number of Seeds Remaining"
  ) +
  my_theme