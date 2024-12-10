
# Fit the age acceleration ~ birth year model with a couple of variations to 
# check whether the random effects structure we used for the model affected our
# results. The data and code are the same as in the paper.

#load packages
library(tidyverse)
library(tidybayes)
library(brms)

# Load aging data
epi_dat <- readRDS('input/PB_clock_ages.rds') %>%
  rename('Year' = yr) %>%
  mutate(Born = Year - round(Age),
         across(c(Year, AgeAccel, Age, Born), 
                list(sc = function(x) as.vector(scale(x, center = T)))))

# Model epigenetic acceleration ~ birth year (model used in paper)
accel_born_mod <-  brm(AgeAccel_sc ~ Born_sc + Sex + (Born_sc + Sex | BearID),
                       data = epi_dat, family = gaussian, 
                       iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                       prior = prior(normal(0,1), class = b),
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       backend = 'cmdstanr')

# Model epigenetic acceleration ~ birth year (model suggested by reviewer)
accel_born_mod2 <-  brm(AgeAccel_sc ~ Born_sc + Sex + (1 | BearID),
                       data = epi_dat, family = gaussian, 
                       iter = 10000, warmup = 5000, chains = 4, cores = 4, 
                       prior = prior(normal(0,1), class = b),
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       backend = 'cmdstanr')
