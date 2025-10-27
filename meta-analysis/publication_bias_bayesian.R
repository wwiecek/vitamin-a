library(rstan)
library(tidybayes)
rstan_options(auto_write = TRUE)

mod <- stan_model("meta-analysis/selection_rubin.stan")
mod2 <- stan_model("meta-analysis/selection_ak.stan")
source("meta-analysis/prepare_ma_data.R")
imdad <- imdad2022 %>%  #for convenience
  filter(group != "Lin 2008")

# Testing Bayesian model on simulated datasets -----

sim <- replicate(100, {
  df <- ma_generator(rel_pp = 0.2) 
  
  ak_result <- try(metastudies_estimation(
    df$yi,
    df$sei,
    cutoffs = 1.96,
    symmetric = TRUE,
    model = "normal"
  ))
  
  if(inherits(ak_result, "try-error"))
    ak_res <- matrix(NA, 2, 3)
  else
    ak_res <- ak_result %>% bind_cols() %>% t() %>% round(3)
  
  fit <- sampling(mod, data = list(y = df$yi, se = df$sei, K = nrow(df), c = 1.96), refresh = 0)
  pars <- c("tau", "sigma_tau", "omega")
  x <- extract(fit, pars)
  bayes_res <- rbind(mean = sapply(x, mean), sd = sapply(x, sd))
  
  rbind(ak_res, bayes_res) %>% cbind(ak = c(1,1,0,0))
})

validrows <- !is.na(sim[1,3,]) & (sim[1,3,] < 10)
sum(validrows)
# summary of point estimates mu, tau, Pr AK parameters:
apply(sim[1,,validrows], 1, summary)
# summary of mean Bayesian parameters
apply(sim[3,,], 1, summary)

# Distribution of omegas
v1 <- sim[1,3,validrows]; v2 <- sim[3,3,]
tibble(val = c(v1, v2), grp = rep(c("AK omega", "Bayes omega"), c(length(v1), length(v2)))) %>%
  ggplot(aes(x = val, fill = grp)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.4, position = "identity", binwidth = 0.2) +
  geom_density(alpha = 0.7) +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 2))



# Fitting to Imdad et al data -----

fit <- sampling(mod, 
                data = compose_data(transmute(imdad, y = tau, se), 
                                    K = nrow(imdad), 
                                    c = 1.96),
                refresh = 0)

fit2 <- baggr(imdad)

print(fit, c("tau", "sigma_tau", "omega"))
fit2

fit3 <- sampling(mod2, refresh = 0,
                 data = compose_data(transmute(imdad, y = tau, se), N = nrow(imdad), c = 1.96))
fit3


# Fun test of water m-a datasets:
load("~/github/water-ma/data/final/ma_datasets.Rdata")
df_main_ma_adj
fitw <- sampling(mod2, refresh = 0,
                 data = compose_data(transmute(df_main_ma_adj, y = tau_f, se = se_f), 
                                     N = nrow(df_main_ma_adj), c = 1.96))
fitw
metastudies_estimation(
  df_main_ma_adj$tau_f,
  df_main_ma_adj$se_f,
  cutoffs = 1.96,
  symmetric = TRUE,
  model = "normal"
) %>% bind_cols() %>% t() %>% round(3)
as.matrix(fitw)[,"omega"] %>% density %>% plot 
