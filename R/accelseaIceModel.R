
# 01 - Prep workspace ====

# Load packages
library(tidyverse)
library(sf)

# Load age data
age_dat <- readRDS('output/new_clock_ages.rds') %>%
  rename(BearID = id) %>%
  mutate(SampleYear = as.numeric(substr(YMD, 1, 4))) %>%
  # For each sample, make a list of the years between the bear's year of
  # birth and its date of sample (these are the days it was alive before the
  # sample)
  rowwise() %>%
  mutate(LiveYears = list(seq(from = Born, to = SampleYear, by = 1)))

# Load ice season data
season_dat <- readRDS('input/NSIDC_ice_season_subpops.rds') %>%
  rename(Population = PopID)


# Combine age and ice-free season data
age_season_dat <- age_dat %>%
  filter(Born  >= 1979) %>%
  select(BearID, sampleId, Population, ageAccel, LiveYears) %>%
  left_join(season_dat, relationship = 'many-to-many') %>%
  filter(Year %in% LiveYears) %>%
  group_by(sampleId) %>%
  # Summarize mean length of ice free seasons from the bear's birth year to its
  # sample date
  summarize(BearID = unique(BearID),
            iceFreeLength = mean(iceFreeLength), ageAccel = unique(ageAccel),
            Population = unique(Population))

# Mean age acceleration by Bear ID (to avoid multiple samples per individual
# and fit simple linear model)
age_season_means <- age_season_dat %>%
  group_by(BearID) %>%
  summarize(iceFreeLength = mean(iceFreeLength), ageAccel = mean(ageAccel),
            Population = unique(Population))

# Plot age by ice-free season length
  age_season_means %>%
    ggplot(aes(x = iceFreeLength, y = ageAccel)) +
    geom_point(aes(colour = Population)) +
    geom_smooth(method = 'lm', colour = 'black', fill = 'black') +
    theme(panel.background = 
            element_rect(fill = 'white', colour = 'black', linewidth = 1),
          plot.margin = 
            unit(c(0.5, 0.5, 1, 1), 'cm'),
          axis.title.x = 
            element_text(size = 15, colour = 'black', vjust = -5),
          axis.title.y = 
            element_text(size = 15, colour = 'black', vjust = 5),
          axis.text =
            element_text(size = 15, colour = 'black'),
          legend.title = element_text(size = 15, colour = 'black'),
          legend.text = element_text(size = 15, colour = 'black'),
          legend.key = element_rect(linewidth = 0.5)
    ) +
    scale_color_viridis_d() +
    labs(x = 'Mean length of ice-free season during lifetime (days)',
         y = 'Epigenetic age acceleration')
  
# Fit linear model
summary(lm(ageAccel ~ iceFreeLength, data = age_season_means))


