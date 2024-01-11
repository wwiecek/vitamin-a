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
regtest(metafor_fe)
regtest(metafor_re)
