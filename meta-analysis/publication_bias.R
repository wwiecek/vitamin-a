# Publication bias assessment for Imdad et al data
library(tidyverse)
library(baggr)
library(metafor)
set.seed(1990)
source("meta-analysis/prepare_ma_data.R")
imdad <- imdad2022 %>%  #for convenience
  filter(group != "Lin 2008")
imdad_nd <- filter(imdad, group != "DEVTA 2013")

metafor_re <- rma(data = imdad, yi = tau, sei = se)
metafor_fe <- rma(data = imdad, yi = tau, sei = se, method = "FE")
funnel(metafor_fe)
funnel(metafor_re)

# Egger's test for funnel plot asymmetry tends to reject null
regtest(metafor_fe)
regtest(metafor_re)
regtest(metafor_fe, model = "lm")
regtest(metafor_fe, model = "lm", predictor = "vi")
# cannot perform Peters' test because I don't have raw numbers of events,
# but unlikely it would change our mind

# Begg and Mazdumar's rank test shows nothing
# (correlation of effect sizes and sampling variances)
# but I think this is a less powerful test
ranktest(metafor_fe)

# Andrews & Kasy publication bias estimate:
# https://maxkasy.github.io/home/metastudy/
# imdad %>% select(tau, se) %>% write_csv("kasy_shiny.csv")
# estimates nonsensical values
# but for symmetrical publication Pr (fewer paramters)
# mu = -0.131 (SE = 0.074), tau = 0.085 (SE = 0.065), Pr = 0.17
# so the penalised result is similar to FE but with predictably larger uncertainty
# (but also appropriately smaller heterogeneity)



# Assessing power of tests -----

# assume we have the same number of studies and the SEs are drawn from empirical set
# true effects follow the RE model we estimated 

