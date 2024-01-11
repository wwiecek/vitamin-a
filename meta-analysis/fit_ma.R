# Re-running meta-analysis models

rm(list = ls())
library(tidyverse)
library(baggr)
library(metafor)
set.seed(1990)
source("meta-analysis/prepare_ma_data.R")
imdad <- imdad2022 %>%  #for convenience
  filter(group != "Lin 2008")
imdad_nd <- filter(imdad, group != "DEVTA 2013")



# Results from metafor -----
rma.uni(yi = tau, sei = se, data = imdad)



# Meta-analysis with baggr -----
bg_imdad_full <- baggr(imdad, pooling = "full", iter = 1e04)
bg_imdad_partial <- baggr(imdad, iter = 1e04)
# No DEVTA:
bg_imdadnd_full <- baggr(imdad_nd, pooling = "full", iter = 1e04)
bg_imdadnd_partial <- baggr(imdad_nd, iter = 1e04)
# No pool:
bg_imdad_nopool <- baggr(imdad, pooling = "none", iter = 1e04)



# Cross-validation models -----

# LOO CV
cv1 <- loocv(imdad, pooling = "full", return_models = T)
cv2 <- loocv(imdad, pooling = "partial", iter = 5000, return_models = T) 
#there should be only a few DTs in each partial pool model but it should be OK
#please check, though, when running this

# LOO CV
cv1n <- loocv(imdad_nd, pooling = "full", return_models = T)
cv2n <- loocv(imdad_nd, pooling = "partial", iter = 5000, return_models = T) 

# Comparison of CV values
data.frame(
  model = 1:nrow(imdad_nd),
  full = cv1n$pointwise,
  partial = cv2n$pointwise) %>% 
  gather(key, value, -model) %>% 
  ggplot(aes(x = model, y = value, color = key, group = key)) + geom_line()
data.frame(
  model = 1:nrow(imdad),
  full = cv1$pointwise,
  partial = cv2$pointwise) %>% 
  gather(key, value, -model) %>% 
  ggplot(aes(x = model, y = value, color = key, group = key)) + geom_line()



# Understanding differences between I^2 and average pooling metric -----

i2 <- pl <- list()
df <- expand.grid(sigma_tau = c(.1, 1, 5),
                  se_i = c(.1, 1, 5),
                  K = c(5, 25, 100))
for(i in 1:nrow(df)) {
  y_i <- rnorm(df$K[i], 0, df$sigma_tau[i])
  y_hat_i <- rnorm(rep(1, df$K[i]), y_i, df$se_i[i])
  dt <- data.frame(tau = y_hat_i, se = df$se_i[i])
  bg <- baggr(dt)
  mf <- rma.uni(yi = tau, sei = se, data = dt)
  pl[[i]] <- pooling(bg)
  i2[[i]] <- summary(mf)$I2
}


# Saving outputs -----

# Save summary statistics (for version control)
bgc <- baggr_compare(
  "Full" = bg_imdad_full,
  "Partial" = bg_imdad_partial,
  "Full, no DEVTA" = bg_imdadnd_full,
  "Partial, no DEVTA" = bg_imdadnd_partial)



write_csv(
  rbind(
    bgc$mean_trt %>% 
      round(2) %>% 
      as.data.frame() %>% 
      rownames_to_column() %>% 
      mutate(measure = "mean"),
    bgc$sd_trt %>% 
      round(2)%>% 
      as.data.frame() %>% 
      rownames_to_column() %>% 
      mutate(measure = "sd"),
    bgc$posteriorpd_trt %>% 
      round(2)%>% 
      as.data.frame() %>% 
      rownames_to_column() %>% 
      mutate(measure = "ppc")
  ),
  "meta-analysis/main_results_check.csv")

# save.image("c:/github/vitamin-a/rerun.Rdata")
save.image("meta-analysis/fit_ma.Rdata")
