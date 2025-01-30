
library(tidyverse)
library(glmnet)
library(brms)

# 1 - Load data and source functions ====

# Load related bears
sibs <- readRDS('input/full_sibs.rds')

# Normalized betas from the original WH population (with sites removed that don't align to 
# the genome and low-quality samples removed)
meth_dat <- readRDS('input/samples_original.rds')

# Get bears with multiple samples
multi_samps <- meth_dat %>%
  group_by(BearID) %>%
  summarize(N = n()) %>%
  filter(N >= 2) %>%
  pull(BearID)

# 3 Run 500 iterations ====
# In each iteration, a test and training sample are taken, the clock is fit and
# validated on the test set, then the test set is combined with additional data
# to test the relationship between birth year and age acceleration

# Set iterations
it <- 1
# Run while loop cross-validation
while(it <= 500) {
  
  print(it)
  
  # Specify training and testing data
  meth_betas_train <- meth_dat %>%
    group_by(Age, Spec, Sex) %>%
    sample_n(1) %>%
    filter(! BearID %in% sibs & ! BearID %in% multi_samps)
  
  meth_betas_test <- meth_dat %>%
    filter(! SampleID %in% meth_betas_train$SampleID)
  
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
  
  # Get betas and ages
  betasLoop <- meth_betas_train_m
  ageLoop <- meth_betas_train$Age
  
  # Make sure ages for training match training samples
  meth_betas_train$chip.ID.loc == rownames(meth_betas_train_m)
  
  # 6 Fit clock and predict on training data ====
  
  # Glmnet model (training betas ~ ages)
  cvfit <- cv.glmnet(betasLoop, ageLoop, nfolds = 10, alpha = .5)
  
  # Add predictions as column to ages in testing data
  preds_test <- meth_betas_test %>%
    select(SampleID:Spec) %>%
    # Predict ages
    mutate(AgePredict = as.numeric(predict(cvfit, newx = meth_betas_test_m, 
                                           type = "response", s = "lambda.min")))
  
  # MAE for testing data (validation)
  MAE <- median(abs(preds_test$AgePredict - preds_test$Age))
  # Pearson's correlation for testing data (validation)
  corr <- as.numeric(cor.test(preds_test$AgePredict, preds_test$Age)$estimate)
  
  
  # Add residuals
  age_accels <- preds_test %>%
    mutate(AgeAccel = lm(preds_test$AgePredict ~ preds_test$Age)$residuals,
           # Add year
           yr = as.numeric(substr(SampleID, 8, 11)),
           Born = yr - floor(Age)) %>%
    # Get mean age acceleration
    group_by(BearID) %>%
    summarize(Born = mean(Born), AgeAccel = mean(AgeAccel))
  
  # Fit model age acceleration ~ year of birth
  m <-  brm(AgeAccel ~ Born,
            data = age_accels, family = gaussian, 
            iter = 10000, warmup = 5000, chains = 4, cores = 4, 
            prior = prior(normal(0,1), class = b),
            control = list(adapt_delta = 0.99, max_treedepth = 20),
            backend = 'cmdstanr')
  
  # Predict across range of birth years (1965-2021)
  # New data
  nd <- data.frame(Born = 1965:2018)
  fitted_m <- fitted(m, probs = c(0.025, 0.975), newdata = nd, summary = F) %>%
    # Make data frame and pivot
    data.frame() %>%
    pivot_longer(everything()) %>%
    bind_cols(expand_grid(draws = 1:20000, nd)) %>%
    group_by(Born) %>%
    # Summarize mean age accel ~ birth year
    summarize(AgeAccel = mean(value)) %>%
    mutate(iteration = it) %>% relocate(iteration)
  
  # Pull out posterior draws
  b <- as.data.frame(m) %>%
    pull(b_Born)
  
  # Add metrics and posterior draws to growing objects
  if(it == 1) {
    # Posterior draws
    b_draws <- b
    # Posterior predictive mean age accel ~ birth year relationship
    b_preds <- fitted_m
    # Clock metrics
    mets <- data.frame(iteration = 1, MAE, corr, N = nrow(meth_betas_train))
  } else {
    b_draws <- c(b_draws, b)
    b_preds <- bind_rows(b_preds, fitted_m)
    mets <- mets %>%
      bind_rows(data.frame(iteration = it, MAE, corr, N = nrow(meth_betas_train)))
  }
  
  # Set new iterations
  it <- it + 1
  
}

# 4 Save objects ====

saveRDS(b_draws, 'output/bootstrap_age_birth_draws.rds')
saveRDS(b_preds, 'output/bootstrap_posterior_pred.rds')
saveRDS(mets, 'output/bootstrap_clock_metrics.rds')

# 5 Summarize the bootstrap clock metrics ====

# Mean MAE and confidence intervals
mean(mets$MAE)
quantile(mets$MAE, c(0.025, 0.975))

# Mean Pearson's correlation and confidence intervals
mean(mets$corr)
quantile(mets$corr, c(0.025, 0.975))

# Mean and credible intervals of age acceleration ~ birth year slope
mean(b_draws)
quantile(b_draws, c(0.025, 0.975))

# 6 Plot the posterior predictions ====

b_preds %>%
  ggplot(aes(x = Born, y = AgeAccel, group = iteration)) +
  geom_line(colour = '#425d9c', alpha = 0.1) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1.25),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -3),
        plot.margin = unit(c(0.25, 0.25, 0.75, 0.75), 'cm'),
        panel.grid = element_line(linewidth = 0.5, colour = '#e5e5e5')) +
  xlab('Year of birth') + ylab('Age acceleration (years)')
