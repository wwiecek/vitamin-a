# Script to print out all the relevant text outputs
library(tidyverse)
library(baggr)
load("fit_ma.Rdata")

print(bg_imdadnd_full, exponent = TRUE)
print(bg_imdad_full, exponent = TRUE)
print(bg_imdadnd_partial, exponent = TRUE)
print(bg_imdad_partial, exponent = TRUE)
print(bg_imdad_nopool, exponent = TRUE)
print(bg_imdad_partial, exponent = FALSE)
heterogeneity(bg_imdad_partial, metric = "i")
heterogeneity(bg_imdadnd_partial, metric = "i")

effect_draw(bg_imdad_partial, trans=exp, s=T) %>% round(2)
sum(effect_draw(bg_imdad_partial, trans=exp) < 1)/20000

cv1
cv2
cv1n
cv2n
loo_compare(cv1, cv2)
loo_compare(cv1n, cv2n)
