load("meta-analysis/fit_ma.Rdata")
library(tidyverse)
library(baggr)
library(metafor)
library(ggdist)
theme_set(theme_minimal(base_size = 10))

# There is practically no impact of prior, i.e.  estimates
# updated with default priors are practically the same as inputs:
forest_plot(bg_imdad_nopool, show = "both")
forest_plot(bg_imdad_partial, show = "both")



# Figure: results of RE model -----
baggr_plot(bg_imdad_partial, hyper = TRUE, transform = exp) + 
  xlab("risk ratio") +
  geom_vline(xintercept = 1, lty = "dashed")
ggsave(file = "figures/baggr-re.pdf", width = 12, height = 12, units = "cm")

baggr_compare("Partial, all data" = bg_imdad_partial, 
              "Full, all data" = bg_imdad_full,
              "Partial, no DEVTA" = bg_imdadnd_partial,
              "Full, no DEVTA" = bg_imdadnd_full,
              "No pooling" = bg_imdad_nopool)

bgc <- baggr_compare("Partial, all data" = bg_imdad_partial, 
                     "Partial, no DEVTA" = bg_imdadnd_partial,
                     compare = "effects") 
bgc
plot(bgc) + ggtitle("Possible treatment effect (new implementation)",
                    "Partial pooling only")

bgc_f <- baggr_compare("Full, all data" = bg_imdad_full,
                       "Full, no DEVTA" = bg_imdadnd_full,
                       compare = "effects") 
bgc_f
plot(bgc_f) + ggtitle("Possible treatment effect (new implementation)",
                      "Full pooling only")



# p.p.c. comparison of RE and FE -----
library(patchwork)
plot(bgc) + 
  ggtitle("Posterior distribution for possible treatment effect", "Partial pooling") +
  plot(bgc_f) + ggtitle("", "Full pooling") +
  plot_layout() &
  theme(legend.position='bottom') &
  scale_x_continuous(limits = c(-2, 1)) &
  xlab("log(RR)")
ggsave(file = "figures/baggr-density.pdf", width = 16, height = 9, units = "cm")



# Figure: out-of-sample p.p.d. (cross-validation models) -----

dt1 <- 
  lapply(cv1$models, effect_draw) %>% 
  as.data.frame() %>% 
  setNames(imdad$group) %>% 
  gather() %>%
  group_by(key) %>% 
  mean_qi()

dt2 <- 
  lapply(cv2$models, effect_draw) %>% 
  as.data.frame() %>% 
  setNames(imdad$group) %>% 
  gather() %>%
  group_by(key) %>% 
  mean_qi()

imdad_ci <- imdad %>% 
  mutate(value = tau,
         key = group,
         .lower = tau - 1.96*se,
         .upper = tau + 1.96*se)

rbind(dt1 %>% mutate(model = "Fixed-effects model (p.p.d.)"), dt2 %>% mutate(model = "Random-effects model (p.p.d.)")) %>% 
  ggplot(aes(y = key, x = value, xmin = .lower, xmax = .upper, group = model, color = model)) + 
  geom_pointinterval(position = position_dodge()) +
  # geom_point(aes(y = group, x = tau), data = mutate(imdad, model = "data"), color = "red")
  geom_pointinterval(data = mutate(imdad_ci, model = "Input data")) +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(-1.5, 1)) +
  xlab("log(RR)") + ylab("") +
  ggtitle("Impact of excluding individual studies on FE and RE estimates") +
  scale_color_discrete(name = "")

ggsave(file = "figures/baggr-oos-ppc.pdf", width = 16, height = 18, units = "cm")
