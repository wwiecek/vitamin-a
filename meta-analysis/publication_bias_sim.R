# Publication bias assessment for Imdad et al data
library(tidyverse)
library(baggr)
library(metafor)
set.seed(1990)
source("meta-analysis/prepare_ma_data.R")

imdad <- imdad2022 %>%  #for convenience
  filter(group != "Lin 2008")
imdad_nd <- filter(imdad, group != "DEVTA 2013")


# Assessing performance of tests -----

# assume we have the same number of studies and the SEs are drawn from empirical set
# true effects follow the RE model we estimated -0.29 (posterior mean)
# with tau = 0.25 (posterior mean)
# relative publication probability is 0.17

ma_generator <- function(se_sim = imdad$se,
                         rel_pp = 0.17,
                         mu = -0.29,
                         tau = 0.25,
                         Ns = 18){
  yi <- vector()
  sei <- vector()
  
  # run studies until we have 18 of them published
  k <- 0
  while(k < Ns){
    rr <- rnorm(1, mu, tau)
    # se <- sample(se_sim, 1)
    se <- runif(1, 0.1, 1)
    rr_hat <- rnorm(1, rr, se)
    z <- abs(rr_hat)/se
    published_pr <- ifelse(z > 1.96, 1, rel_pp)
    pub <- rbinom(1, 1, published_pr)
    if(pub) {
      k <- k+1
      yi[k]  <- rr_hat
      sei[k] <- se
    }
  }
  data.frame(yi, sei)
}

# Now compare a few different settings (by hand!) for ma_generator():

runsim <- function(...)
  replicate(100, {
    df <- ma_generator(...) 
    ak_result <- try(metastudies_estimation(
      df$yi, 
      df$sei, 
      cutoffs = 1.96, 
      symmetric = TRUE, 
      model = "normal"
    ))
    if(inherits(ak_result, "try-error"))
      return(matrix(NA, 2, 3))
    ak_result %>%
      bind_cols() %>%
      t() %>% 
      round(3)
  })

# summary of mu, tau, Pr:
sim1 <- runsim(rel_pp = .17)
sim2 <- runsim(rel_pp = 1)


apply(sim[1,,], 1, summary)
sim[1,3,] %>% hist(breaks = 100, xlim = c(0, 2))


# This is what I saw:
# it performs quite nicely given that it's only 18 studies
# i.e. it retrieves the parameters quite well if the model is correctly specified

# if pr of publication is 1 it often fails and the probability is just extremely unstable
# because the method cannot estimate it

# Egger's test:
sim_egg <- replicate(1000, {
  df <- ma_generator(mu = -0.131, tau = 0.085) 
  rm <- try(rma(yi=yi, sei=sei, data = df))
  if(inherits(rm, "try-error"))
    return(NA)
  rt <- regtest(rm)
  rt$pval
})
hist(sim_egg)
sum(sim_egg < 0.05)/1000
# low power to find funnel plot asymmetry under this model
# (because this model does not work on discontinuity)
