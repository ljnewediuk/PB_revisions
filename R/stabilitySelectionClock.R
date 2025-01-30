# Stability selection clock

# Use stability selection to find sites consistently selected
# Re-do this using training and testing data, then use the stably selected CpGs
# to predict age in testing data vs. a classic LASSO model

#  1 - Prep workspace ====

# Load libraries
library(tidyverse)
library(glmnet)

# Load functions needed for stability selection
source('functions/stabSelFunctions.R')

# 2 - Load data ====

# Normalized betas from all populations (with sites removed that don't align to 
# the genome and low-quality samples removed)
# This file is available through MS OneDrive as it's too large to push to GitHub.
meth_dat <- readRDS('input/samples_new.rds')

# Names of cpg sites
cpg_names <- meth_dat %>%
  select(! sampleId:Population) %>% 
  colnames()

# 3 - Prep data for clock ====

# Check which bears have repeat samples
multi_samps <- meth_dat %>%
  group_by(id) %>%
  summarize(n_samps = n()) %>%
  filter(n_samps > 2) %>%
  pull(id)

# Set seed
set.seed(1234)

# Separate into training, testing, validation sets 
# (sample by population, tissue, age)
train_dat <- meth_dat %>%
  filter(! id %in% c(multi_samps)) %>%
  group_by(Population, Spec, age) %>%
  slice_sample(n = 1)

# Testing data distinct from training data
test_dat <- meth_dat %>%
  filter(! sampleId %in% train_dat$sampleId) %>%
  distinct()

# Get ages from data for feature selection (vector length = nrow(DNAm))
fs_ages <- as.numeric(train_dat$age)

# Get matrix of DNAm data
fs_meth <- train_dat %>%
  # Make chip positions rownames
  column_to_rownames('chip.ID.loc') %>%
  # Remove extra cols
  select(! sampleId:Population) %>%
  # Convert to matrix
  as.matrix()

# Stability selection

# Fit LASSO clock
mod.cv <- cv.glmnet(x = fs_meth, y = fs_ages, alpha = 1, nfolds = 10)

# CpGs from lasso clock (39 sites)
lasso_clock <- as.matrix(coef(mod.cv, s = 'lambda.min')) %>%
  as.data.frame() %>%
  rownames_to_column('cg') %>%
  dplyr::rename('beta' = s1) %>%
  filter(beta != 0)

# Run stability selection x1000
result <- lapply(1:1000, FUN = stabsel)

# Compute average number of selected variables (q_stab)
p_list <- list()
for ( i in 1 :length(result)) {
  # Get length of each result (number of selected CpGs)
  res <- as.numeric(result[[i]])
  p_list[[i]] <- length(res)
}
p_out <- do.call(cbind,p_list)
# Average number of selected CpGs
q_stab <- mean(p_out)

# Compute the selection probabilities

n_select <- rep(0, length(cpg_names))
t_select <- 0 * n_select

for (i in 1:length(result)) {
  # Get names of CpGs in result i
  tt <- names(result[[i]])[-1]
  # Add one for each selected CpG
  if(length(tt)==0){
    
  } else {
    t_select[which(cpg_names %in% names(result[[i]])[-1])] <- 1
    n_select <- n_select + t_select
    # Re-zero
    t_select <- 0*n_select
  }
}

# Get proportion of times each site was selected
out <- n_select/length(result)
df_out <- data.frame(cpg_names = cpg_names, pi_select = out)
# Arrange in descending order
select_prob <- df_out %>%
  arrange(desc(pi_select))

# Find selection probability threshold (max 2 false discoveries; threshold = E_v)
thresh <- finding_thresh(q = q_stab, p = length(cpg_names), E_v = 2)

# Number of stably predictive CpGs (5)
n_cpgs <- nrow(select_prob[which(select_prob$pi_select > thresh),])

# Fit clock with top 5 stably-selected CpGs
i <- 5
mod_best <- ssGAM(train_dat, test_dat, select_prob, addPlot = T)

# Predict ages
preds <- test_dat %>%
  mutate(agePred = mod_best$data$agePredict) %>%
  select(sampleId:Population, agePred)

# Add residuals (age acceleration)
preds$ageAccel <- lm(agePred ~ age, data = preds)$resid

# Plot the clock
preds %>%
  ggplot(aes(x = age, y = agePred, colour = Spec)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_point() +
  geom_smooth(method = 'lm', colour = 'black') +
  ylim(-1, 30) + xlim(-1, 30) +
  scale_colour_manual(values = c('#c2c973', '#da7959', '#8ba5a9')) +
  theme_bw() +
  labs(y = 'Epigenetic age', x = 'Chronological age')

# Median absolute error
mod_best$mae
# Correlation
mod_best$cor
# R^2
mod_best$R2

# Filter out just the western Hudson Bay population (used in original paper)
WH_preds <- preds %>% 
  filter(Population == 'WH') %>%
  # Mean age acceleration by individual (to deal with the few individuals having
  # multiple samples each)
  group_by(id) %>%
  summarize(Born = unique(Born), ageAccel = mean(ageAccel))

# Fit linear model
summary(lm(ageAccel ~ Born, data = WH_preds))

# Plot
WH_preds %>%
  ggplot(aes(x = Born, y = ageAccel)) + 
  geom_point() + 
  geom_smooth(method = 'lm', colour = 'black') +
  theme_bw() +
  labs(y = 'Epigenetic age acceleration (residuals)', x = 'Birth year')

# Save the age predictions
saveRDS(preds, 'output/new_clock_ages.rds')
