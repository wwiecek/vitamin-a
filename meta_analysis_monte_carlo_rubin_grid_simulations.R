# Script to generate Meta analysis Monte Carlos Simulations
# Author: Hannah Balikci, based on code by Rachael Meager
# Date: August 2023

# Running meta_analysis_monte_carlo_rubin_grid_simulations.R 
# will save the relevant outputs in simulations_results/

### PRELIMINARIES ###

rm(list = ls())

library(multiwayvcov)
library(sandwich)
library(stargazer)
library(foreign)
library(Hmisc)
library(xtable)
library(coda)
library(devtools)
library(gridExtra)
library(ggplot2)
library(parallel)
library(rstan)

# function to convert stan to coda
stan2coda <- function(fit2) {
  mcmc.list(lapply(1:ncol(fit2), function(x) mcmc(as.array(fit2)[,x,])))
}

## GENERATE THE DATA ###

# set the seed
set.seed(20)

# meta-parameters and storage
S <- 5 # number of simulations
K <- 9 # number of studies

bhm_tau_error <- matrix(rep(NA,16),nrow=4,ncol=4)
fe_tau_error <- matrix(rep(NA,16),nrow=4,ncol=4)

# these define the grid cells for the comparative performance evaluations

se_min <- c(0.1,5,10,15)
se_max <- c(5,10,15,20)
sigma_min <- c(0.1,5,10,15)
sigma_max <- c(5,10,15,20)

# Call the stan model "rubin_model_code.stan" which must be placed in same directory

model_file <- file.path('/Users/hbalikci/Documents/Chicago_Harris/Life Admin/DIL/Vitamin A/vitamin-A_dropbox/rubin_model_code.stan')
model <- stan_model(model_file)

sim_types <- c("student_t_16_cells", "location_outlier_16_cells", "precision_outlier_16_cells", "normal")

# Loop through simulation types
for (sim_type in sim_types) {
  for(i in 1:4){
    for(j in 1:4){
      
      tau_k_draws <- matrix(rep(0,K*S), nrow=S, ncol=K) # k-specific effects
      hat_tau_k_draws <- matrix(rep(0,K*S), nrow=S, ncol=K) # k-specific estimates
      se_k_draws <- matrix(rep(0,K*S), nrow=S, ncol=K) #k-specific standard errors
      
      seed <- (i+1)*(j+1)
      
      tau_draws <- rnorm(n=S,mean=3,sd=1) # parent mean!
      sigma_tau_draws <- runif(n=S,min=sigma_min[i],max=sigma_max[i]) # parent standard deviation!
      
      for (s in 1:S) {
        seed <- 2*(i+1)*(j+1)
        
        if (sim_type == "student_t_16_cells") {
          tau_k_draws[s,] <- tau_draws[s] + sqrt(sigma_tau_draws[s])*rt(n=K, df=3)
          se_k_draws[s,] <- runif(n=K, min=se_min[j], max = se_max[j])
        } else if (sim_type == "location_outlier_16_cells") {
          tau_k_draws[s,] <- rnorm(n=K, mean=tau_draws[s], sd=sigma_tau_draws[s])
          tau_k_draws[s, 1] <- rnorm(n=1, mean=tau_draws[s]+10, sd=sigma_tau_draws[s])
          se_k_draws[s,] <- runif(n=K, min=se_min[j], max=se_max[j])
        } else if (sim_type == "precision_outlier_16_cells") {
          tau_k_draws[s,] <- rnorm(n=K, mean=tau_draws[s], sd=sigma_tau_draws[s])
          se_k_draws[s,] <- runif(n = K, min = se_min[j], max = se_max[j])
          se_k_draws[s, 1] <- runif(n=1, min=0.1*se_min[j], max=0.1*se_max[j])
        } else if (sim_type == "normal") {
          tau_k_draws[s,] <- rnorm(n=K, mean=tau_draws[s], sd=sigma_tau_draws[s])
          se_k_draws[s,] <- runif(n=K, min=se_min[j], max=se_max[j])
        }
        
        for (k in 1:K) {
          hat_tau_k_draws[s, k] <- rnorm(n = 1, mean = tau_k_draws[s, k], sd = se_k_draws[s, k])
        }
      }
      
      ## RUN THE ANALYSIS ##
      
      # FE-analysis
      
      # define function
      FE_estimator <- function(taus, ses){
        inverse_vars <- 1/ (ses^2)
        inverse_var_tot <- sum(inverse_vars)
        inverse_var_weight <- inverse_vars/inverse_var_tot
        FE_mean <- sum(taus* inverse_var_weight)
        FE_se <- sqrt(1/inverse_var_tot)
        FE_output <- c(FE_mean,FE_se)
        return(FE_output)
      }
      
      # storage
      FE_means <- vector(mode="numeric",length=S)
      FE_ses <- vector(mode="numeric",length=S)
      
      for (s in 1:S){
        
        FE_means[s] <- FE_estimator(hat_tau_k_draws[s,],se_k_draws[s,])[1]
        FE_ses[s] <-  FE_estimator(hat_tau_k_draws[s,],se_k_draws[s,])[2]
        
      } #closes the for-loop indexed by s
      
      
      # evaluate performance of means 
      FE_means_error <- FE_means - tau_draws
      FE_means_mse <- mean(FE_means_error^2)
      
      fe_tau_error[i,j] <- FE_means_mse
      print(fe_tau_error[i,j])
      
      ### BAYESIAN HIERARCHICAL ANALYSIS ###
      
      #storage
      stan_output_summaries <- list("eg","eg","eg")
      posterior_mean_tau <- vector(mode="numeric", length=S)
      posterior_sd_tau <- vector(mode="numeric", length=S)
      posterior_sigma_tau <- vector(mode="numeric",length=S)
      posterior_mean_tau_ks <- matrix(rep(0,K*S), nrow=S, ncol=K) # k-specific estimates
      
      
      # this little bit is J Huggins' genius set up #
      # some knobs we can tweak
      chains <- 4
      iters <- 5000
      control <- list(adapt_t0 = 10,       # default = 10
                      stepsize = 1,        # default = 1
                      max_treedepth = 10)   # default = 10
      seed <- 3*(i+1)*(j+1)
      
      # now loop through the data sets and run the model on each
      for(s in 1:S){
        dataset_temp <- list(tau_hat_k = hat_tau_k_draws[s,],
                             se_k = se_k_draws[s,],
                             K = K)
        
        sflist_temp <-
          mclapply(1:chains, mc.cores = chains,
                   function(i) sampling(model, data = dataset_temp, seed = seed,
                                        chains = 1, chain_id = i, # refresh = -1,
                                        iter = iters, control = control))
        stanfit_temp <- sflist2stanfit(sflist_temp)
        
        summary_stanfit_temp <- summary(stanfit_temp)
        textable_stanfit_temp <- xtable(summary_stanfit_temp$summary)
        posterior_mean_tau[s] <- textable_stanfit_temp[1,1]
        posterior_sd_tau[s] <- textable_stanfit_temp[1,3]
        posterior_sigma_tau[s] <- textable_stanfit_temp[2,1]
        posterior_mean_tau_ks[s,] <- textable_stanfit_temp[3:(K+2),1]
        print(textable_stanfit_temp)
        stan_output_summaries[[s]] <- textable_stanfit_temp
        rm(stanfit_temp)
        
      } # closes the forloop indexed by s
      
      # now compute errors
      # evaluate performance of means 
      BHM_means_error <- posterior_mean_tau - tau_draws
      BHM_means_mse <- mean(BHM_means_error^2)
      
      bhm_tau_error[i,j] <- BHM_means_mse
      print(bhm_tau_error[i,j])
      
    } # close the j-loop
  } # close the i-loop
  
  # Print and save results
  bhm_error_table <- data.frame(bhm_tau_error)
  rownames(bhm_error_table) <- table_rownames
  colnames(bhm_error_table) <- table_colnames
  
  fe_error_table <- data.frame(fe_tau_error)
  rownames(fe_error_table) <- table_rownames
  colnames(fe_error_table) <- table_colnames
  
  print(fe_error_table, digits = 3)
  print(bhm_error_table, digits = 3)
  
  # Save objects
  save(fe_error_table, bhm_error_table, file = paste0(sim_type, "_tables.RData"))
  pdf("fe_error_table", sim_type, ".pdf")
  grid.table(data)
  dev.off()
  
  pdf("bhm_error_table", sim_type, ".pdf")
  grid.table(data)
  dev.off()
  
  save.image(file = paste0("meta_analysis_monte_carlo_rubin_grid_if_se_versus_sigma_", sim_type, ".RData"))
} #close sim_type loop

