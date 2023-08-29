# Script to generate meta analysis Monte Carlos Simulations
# Author: Rachael Meager, code reorganised by Hannah Balikci
# Date: August 2023

# This code will require 64 x 5000 models, so you have to be very patient
# It will probably take about a day, depending on a machine?

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
library(gridExtra)

# function to convert stan to coda
stan2coda <- function(fit2) {
  mcmc.list(lapply(1:ncol(fit2), function(x) mcmc(as.array(fit2)[,x,])))
}

## GENERATE THE DATA ###

# set the seed
set.seed(20)

# meta-parameters and storage
S <- 5000 # number of simulations
K <- 9 # number of studies

bhm_tau_error <- matrix(rep(NA,16),nrow=4,ncol=4)
fe_tau_error <- matrix(rep(NA,16),nrow=4,ncol=4)

# Define the grid cells for the comparative performance evaluations
se_min <- c(0.1,5,10,15)
se_max <- c(5,10,15,20)
sigma_min <- c(0.1,5,10,15)
sigma_max <- c(5,10,15,20)

table_rownames <- c("Sigma : 0-5","5-10","10-15","15-20")
table_colnames <- c("SE : 0-5","5-10","10-15","15-20")

# Compile the stan model
rstan_options(auto_write = TRUE)
model <- stan_model('rubin_model_code.stan')

# Define the desired file path for outputs
desired_path_tables <- "figures/"
desired_path <- "simulation_results/"

sim_types <- c("student_t_16_cells", 
               "location_outlier_16_cells", 
               "precision_outlier_16_cells", 
               "normal")

# Loop through simulation types
for (t  in sim_types) {
  for(i in 1:4){
    for(j in 1:4){
      
      tau_k_draws <- matrix(rep(0,K*S), nrow=S, ncol=K) # k-specific effects
      hat_tau_k_draws <- matrix(rep(0,K*S), nrow=S, ncol=K) # k-specific estimates
      se_k_draws <- matrix(rep(0,K*S), nrow=S, ncol=K) # k-specific standard errors
      
      seed <- (i+1)*(j+1)
      
      tau_draws <- rnorm(n=S,mean=3,sd=1) # parent mean
      sigma_tau_draws <- runif(n=S,min=sigma_min[i],max=sigma_max[i]) # parent standard deviation
      
      for (s in 1:S) {
        seed <- 2*(i+1)*(j+1)
        
        if (t == "student_t_16_cells") {
          tau_k_draws[s,] <- tau_draws[s] + sqrt(sigma_tau_draws[s])*rt(n=K, df=3)
          se_k_draws[s,] <- runif(n=K, min=se_min[j], max = se_max[j])
        } else if (t == "location_outlier_16_cells") {
          tau_k_draws[s,] <- rnorm(n=K, mean=tau_draws[s], sd=sigma_tau_draws[s])
          tau_k_draws[s, 1] <- rnorm(n=1, mean=tau_draws[s]+10, sd=sigma_tau_draws[s])
          se_k_draws[s,] <- runif(n=K, min=se_min[j], max=se_max[j])
        } else if (t == "precision_outlier_16_cells") {
          tau_k_draws[s,] <- rnorm(n=K, mean=tau_draws[s], sd=sigma_tau_draws[s])
          se_k_draws[s,] <- runif(n = K, min = se_min[j], max = se_max[j])
          se_k_draws[s, 1] <- runif(n=1, min=0.1*se_min[j], max=0.1*se_max[j])
        } else if (t == "normal") {
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
      # print(fe_tau_error[i,j])
      
      ### BAYESIAN HIERARCHICAL ANALYSIS ###
      
      #storage
      stan_output_summaries <- list("eg","eg","eg")
      posterior_mean_tau <- vector(mode="numeric", length=S)
      posterior_sd_tau <- vector(mode="numeric", length=S)
      posterior_sigma_tau <- vector(mode="numeric",length=S)
      posterior_mean_tau_ks <- matrix(rep(0,K*S), nrow=S, ncol=K) # k-specific estimates
      
      
      # this little bit is J Huggins' genius set up #
      # some knobs we can tweak
      nchains <- 4
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
        
        stanfit_temp <- sampling(model, data = dataset_temp, seed = seed,
                                 chains = nchains, chain_id = i, refresh = 0,
                                 iter = iters, control = control)
        
        summary_stanfit_temp <- summary(stanfit_temp)
        textable_stanfit_temp <- xtable(summary_stanfit_temp$summary)
        posterior_mean_tau[s] <- textable_stanfit_temp[1,1]
        posterior_sd_tau[s] <- textable_stanfit_temp[1,3]
        posterior_sigma_tau[s] <- textable_stanfit_temp[2,1]
        posterior_mean_tau_ks[s,] <- textable_stanfit_temp[3:(K+2),1] # adjusted to allow for higher K
        # print(textable_stanfit_temp)
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
  
  # Save results
  bhm_error_table <- data.frame(bhm_tau_error)
  rownames(bhm_error_table) <- table_rownames
  colnames(bhm_error_table) <- table_colnames
  
  fe_error_table <- data.frame(fe_tau_error)
  rownames(fe_error_table) <- table_rownames
  colnames(fe_error_table) <- table_colnames
  
  # Convert data frames to tables
  fe_table <- tableGrob(format(round(fe_error_table, 3), nsmall = 3))
  bhm_table <- tableGrob(format(round(bhm_error_table, 3), nsmall = 3))
  
  # Define the base filename for the tables
  base_filename <- paste0(t, "_tables")
  
  fe_pdf_path <- file.path(desired_path_tables, paste0(base_filename, "_fe.pdf"))
  ggsave(fe_pdf_path, plot = fe_table, device = "pdf")
  
  bhm_pdf_path <- file.path(desired_path_tables, paste0(base_filename, "_bhm.pdf"))
  ggsave(bhm_pdf_path, plot = bhm_table, device = "pdf")
  
  # Save the simulation data as RData files
  cat("\n", i, j)
  filename <- paste0(t, ".RData")
  full_file_path <- file.path(desired_path, filename)
  save.image(file = full_file_path)
  
} #close sim_type loop
