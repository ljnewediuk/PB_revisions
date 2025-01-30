
# Fit the clock over a series of seeds 1:100 to check stability of the accuracy
# metrics. This is the same code and data we used in the original paper, with
# the addition of the simulation.

# load packages
library(tidyverse)
library(glmnet)

# Training betas
meth_betas_train <- readRDS('input/train_samples_original.rds')

# Testing betas
meth_betas_test <- readRDS('input/test_samples_original.rds')

# List IDs of training and testing bears
train_bears <- meth_betas_train %>%
  pull(BearID)
test_bears <- meth_betas_test %>%
  pull(BearID)

# Get matrix of betas for training data
meth_betas_train_m <- meth_betas_train %>%
  # Make chip positions rownames
  column_to_rownames('chip.ID.loc') %>%
  # Remove extra cols
  select(! SampleID:Spec) %>%
  # Convert to matrix
  as.matrix()

# Get matrix of betas for test data
meth_betas_test_m <- meth_betas_test %>%
  # Make chip positions rownames
  column_to_rownames('chip.ID.loc') %>%
  # Remove extra cols
  select(! SampleID:Spec) %>%
  # Convert to matrix
  as.matrix()

# Check ages

# Add ages for training and testing
age_df <- bind_rows(meth_betas_train, meth_betas_test) %>%
  mutate(AgePredict = 0) %>%
  select(SampleID:chip.ID.loc, AgePredict)

# Get betas and ages
betasLoop <- meth_betas_train_m
ageLoop <- as.numeric(age_df[age_df$SampleID %in% meth_betas_train$SampleID ,]$Age)

# Fit the clocks for random seeds 1:100

mod_df <- data.frame()

for(rseed in 1:100) {
  
  # Set the seed
  set.seed(rseed)
  
  # Glmnet model (training betas ~ ages)
  cvfit <- cv.glmnet(betasLoop, ageLoop, nfolds = 10, alpha = .5)
  
  # Add predictions as column to ages in training data
  age_preds <- age_df %>%
    filter(SampleID %in% meth_betas_test$SampleID) %>%
    # Predict model
    mutate(AgePredict = as.numeric(predict(cvfit, newx = meth_betas_test_m, 
                                           type = "response", s = "lambda.min")))
  # Add residuals
  age_preds <-  age_preds %>%
    mutate(AgeAccel = lm(age_preds$AgePredict ~ age_preds$Age)$residuals,
           # Add year
           yr = as.numeric(substr(SampleID, 8, 11))) %>%
    # Remove chip column
    select(! chip.ID.loc)
  
  # Calculate median absolute error, correlation, R^2, n sites
  age_mae <- median(abs(age_preds$AgePredict - age_preds$Age))
  age_corr <- as.numeric(cor.test(age_preds$AgePredict, age_preds$Age)$estimate)
  age_R2 <- summary(lm(age_preds$AgePredict ~ age_preds$Age))$r.squared
  nSites <- length(coef(cvfit, s = 'lambda.min')[coef(cvfit, s = 'lambda.min')[,1]!= 0])
  
  mod_df_row <- data.frame(seed = rseed, MAE = age_mae, 
                           corr = age_corr, R2 = age_R2, nSites)
  mod_df <- bind_rows(mod_df, mod_df_row)
  
}

# Table of clock accuracy metrics across seeds
mod_df %>%
  pivot_longer(MAE:nSites, names_to = 'metric') %>%
  group_by(metric) %>%
  summarize(mean = mean(value), se = sd(value)/sqrt(100))
