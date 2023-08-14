library(readxl)
library(tidyverse)
library(baggr)
library(metafor)

set.seed(1990)
theme_set(theme_minimal(base_size = 10))
# baggr_theme_update(base_size = NULL)
# load("c:/github/vitamin-a/rerun.Rdata")
load("rerun.Rdata")

# Check intervals in input data -----
mayo <- read_excel("input/mayo_wilson_fig3.xlsx") %>%
  mutate(uci = log(uci) - log(RR), lci = log(lci)- log(RR), RR = log(RR))

# Awasthi numbers are in Fig 4 here:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3647148/


# Note, that for Dibley 1996 there is something going on,
# probably this 0.01 is not precise enough
ggplot(mayo, aes(x = group, y = 0, ymin = lci, ymax=uci)) + 
  geom_point() + geom_errorbar() + coord_flip()


# I think LCI for Dibley study should be 0.0044
exp(-1.11 - 4.3)
# but maybe they didn't want to round to 0.00...


# You can see here how there is .01-.02 error in SE due to rounding off
# of LCI/UCI... but for Dibley it is .07
mayo <- read_excel("input/mayo_wilson_fig3.xlsx") %>%
  mutate(tau = log(RR), 
         se = (log(uci) - log(lci)) / (2*1.96),
         se_alt = (log(uci) - log(RR)) / (1.96))
# We use 'se', not 'se_alt' because they are slightly wider
# so, for me, more trustworthy

mayo_nd <- filter(mayo, group != "DEVTA 2013")



# Results from metafor -----

rma.uni(yi = tau, sei = se, data = mayo)

# Analysis with baggr -----
bg_mayo_full <- baggr(mayo, pooling = "full", iter = 1e04)
bg_mayo_partial <- baggr(mayo, iter = 1e04)
# No DEVTA:
bg_mayond_full <- baggr(mayo_nd, pooling = "full", iter = 1e04)
bg_mayond_partial <- baggr(mayo_nd, iter = 1e04)
# No pool:
bg_mayo_nopool <- baggr(mayo, pooling = "none", iter = 1e04)

# There is practically no impact of prior, i.e.  estimates
# updated with default priors are practically the same as inputs:
forest_plot(bg_mayo_nopool, show = "both")
forest_plot(bg_mayo_partial, show = "both")

print(bg_mayond_full, exponent = TRUE)
print(bg_mayo_full, exponent = TRUE)
print(bg_mayond_partial, exponent = TRUE)
print(bg_mayo_partial, exponent = TRUE)
print(bg_mayo_nopool, exponent = TRUE)

heterogeneity(bg_mayo_partial, metric = "i")
heterogeneity(bg_mayond_partial, metric = "i")

baggr_plot(bg_mayo_partial, hyper = TRUE, transform = exp)
ggsave(file = "figures/baggr-re.pdf", width = 12, height = 12, units = "cm")

baggr_compare("Partial, all data" = bg_mayo_partial, 
              "Full, all data" = bg_mayo_full,
              "Partial, no DEVTA" = bg_mayond_partial,
              "Full, no DEVTA" = bg_mayond_full,
              "No pooling" = bg_mayo_nopool)

bgc <- baggr_compare("Partial, all data" = bg_mayo_partial, 
                     "Partial, no DEVTA" = bg_mayond_partial,
                     compare = "effects") 
bgc
plot(bgc) + ggtitle("Possible treatment effect (new implementation)",
                    "Partial pooling only")

bgc_f <- baggr_compare("Full, all data" = bg_mayo_full,
                       "Full, no DEVTA" = bg_mayond_full,
                       compare = "effects") 
bgc_f
plot(bgc_f) + ggtitle("Possible treatment effect (new implementation)",
                      "Full pooling only")





# Table for the partially pooled model results


library(patchwork)

plot(bgc) + 
  ggtitle("Posterior distribution for possible treatment effect", "Partial pooling") +
  plot(bgc_f) + ggtitle("", "Full pooling") +
  plot_layout() &
  theme(legend.position='bottom') &
  scale_x_continuous(limits = c(-2, 1)) &
  xlab("Treatment effect (log scale)")
ggsave(file = "figures/baggr-density.pdf", width = 16, height = 9, units = "cm")

# LOO CV
cv1 <- loocv(mayo, pooling = "full", return_models = T)
cv2 <- loocv(mayo, pooling = "partial", iter = 5000, return_models = T) 
#there should be only a few DTs in each partial pool model but it should be OK
#please check, though, when running this
loo_compare(cv1, cv2)

# LOO CV
cv1n <- loocv(mayo[-17,], pooling = "full", return_models = T)
cv2n <- loocv(mayo[-17,], pooling = "partial", iter = 5000, return_models = T) 
loo_compare(cv1n, cv2n)
data.frame(
  model = 1:16,
  full = cv1n$pointwise,
  partial = cv2n$pointwise) %>% 
  gather(key, value, -model) %>% 
  ggplot(aes(x = model, y = value, color = key, group = key)) + geom_line()

#there should be only a few DTs in each partial pool model but it should be OK
#please check, though, when running this

# Let's plot each distribution and then the remaining observation:
lapply(cv1n$models, effect_draw)
lapply(cv2n$models, effect_draw)

library(ggdist)

dt1 <- 
  lapply(cv1$models, effect_draw) %>% 
  as.data.frame() %>% 
  setNames(mayo$group) %>% 
  gather() %>%
  group_by(key) %>% 
  mean_qi()

dt2 <- 
  lapply(cv2$models, effect_draw) %>% 
  as.data.frame() %>% 
  setNames(mayo$group) %>% 
  gather() %>%
  group_by(key) %>% 
  mean_qi()

rbind(dt1 %>% mutate(model = "A"), dt2 %>% mutate(model = "B")) %>% 
  ggplot(aes(y = key, x = value, group = model, color = model)) + 
  geom_pointinterval(aes(xmin = .lower, xmax = .upper), position = position_dodge()) +
  geom_point(aes(y = group, x = tau), data = mutate(mayo, model = "data"), color = "red")

# Understand differences between I^2 and average pooling metric -----

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

# save.image("c:/github/vitamin-a/rerun.Rdata")
save.image("rerun.Rdata")
