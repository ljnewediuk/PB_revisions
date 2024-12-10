
# Function for stability selection
stabsel <- function(i) {
  set.seed(i)
  # select subsample (indices)
  ind_sub <- sample (1:length(fs_ages), 
                     size = floor(length(fs_ages)/2), 
                     replace = FALSE)
  # Subset methylation data and ages by index
  fs_meth_sub <- fs_meth[ind_sub,]
  fs_ages_sub <- fs_ages[ind_sub]
  # run LASSO on subsample
  res_sub <- glmnet(x = fs_meth_sub, y = fs_ages_sub, 
                    standardize = T, alpha = 1, lambda = mod.cv$lambda.1se)
  # extract CpGs selected for prediction
  res <- coef(res_sub)
  res <- res[which(res[,1] != 0) ,]
  # Output
  return(res)
}

# Function to compute the threshold (formula from Haftorn et al. 2023)
finding_thresh <- function(q, p, E_v) {
  pi_thresh = ((q^2) / (p * E_v)) + 1 / 2
  return(min(pi_thresh, 1))
}

# Function to fit stability selection gams, plot clocks, and measure accuracy
ssGAM <- function(train_dat, test_dat, prob_dat, thresh = NULL, by_thresh = F, addPlot = F) {
  # Subset results above threshold if by_thresh == T
  if(by_thresh == T) {
    stable_set <- prob_dat[which(prob_dat$pi_select > thresh),]
  } else {
    stable_set <- prob_dat[1:i,]
  }
  # Model
  gam_mod <- mgcv::gam(reformulate(stable_set$cpg_names, response = 'age'), 
                       data = train_dat)
  #  Predict
  age_df <- data.frame(age = test_dat$age, agePredict = predict(gam_mod, test_dat))
  # Calculate mae 
  mae <- as.numeric(ie2misc::mae(age_df$agePredict, age_df$age))
  # And correlation
  cor <- as.numeric(cor.test(age_df$agePredict, age_df$age)$estimate)
  # And R2
  R2 <- summary(lm(agePredict ~ age, data = age_df))$r.squared
  # Number of CpGs
  n <- nrow(stable_set)
  # Plot
  if(addPlot == TRUE) {
    plt <- ggplot(age_df, aes(x = age, y = agePredict)) +
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', linewidth = 1) +
      geom_point(colour = 'black') +
      theme(panel.background = element_rect(fill = 'white',
                                            colour = 'black',
                                            linewidth = 1),
            plot.margin = unit(c(0.5, 0.5, 0.75, 0.75), 'cm'),
            panel.grid = element_blank(),
            axis.title.x = element_text(colour = 'black', size = 15, vjust = -3),
            axis.title.y = element_text(colour = 'black', size = 15, vjust = 5),
            axis.text = element_text(colour = 'black', size = 15),
            legend.position = 'inside',
            legend.position.inside = c(.15,.85),
            legend.title = element_blank(),
            legend.text = element_text(colour = 'black', size = 15),
            legend.key = element_rect(colour = 'black', size = 0.5)) +
      xlim(-2, 30) + ylim(-2, 30) +
      labs(x = 'Chronological age', y = 'Epigenetic age')
  } else {
    plt <- NULL
  }
  # Combine a list with plot, age data, and accuracy measures
  age_list <- list(data = age_df, plot = plt, n_cpgs = n, mae = mae, cor = cor, R2 = R2)
  # Return list
  return(age_list)
}
