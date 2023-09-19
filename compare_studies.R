# Let us compare three different meta-analysis specifications

library(tidyverse)
library(baggr)
library(metafor)
library(ggdist)
theme_set(theme_minimal(base_size = 10))
set.seed(1990)
source("prepare_ma_data.R")

# We find only minor disagreements between the three studies in terms of input data
rbind(
  mayo %>% select(group, tau, se) %>% mutate(source = "Mayo"),
  mayo %>% select(group, tau, se) %>% mutate(source = "Mayo (alternative SE)"),
  awasthi %>% select(group, tau, se) %>% mutate(source = "Awasthi"),
  imdad2022 %>% select(group, tau, se) %>% mutate(source = "Imdad 2022")
) %>% 
  mutate(lci = exp(tau - 1.96*se),
         uci = exp(tau + 1.96*se)) %>% 
  mutate(rr = exp(tau)) %>% 
  ggplot(aes(x = group, y = rr, ymin = lci, ymax=uci, group = source, color = source)) + 
  geom_pointinterval(position = position_dodge()) +
  # geom_point() + 
  # geom_errorbar() + 
  coord_flip() +scale_y_log10() +
  theme(legend.position = "bottom")


# THE REST OF THIS SCRIPT IS UNFINISHED FOR NOW, BUT WE CAN SEE STRONG SIMILARITIES

# How to define standard errors?
ggplot(awasthi, aes(x = group, y = rr, ymin = lci, ymax=uci)) + 
  geom_errorbar(aes(x = group, y = rr, ymin = lci2, ymax=uci2), color="red") +
  geom_point() + 
  geom_point(aes(y=exp(tau)), color="red")+
  geom_errorbar() + coord_flip() +scale_y_log10()

bg_a <- baggr(awasthi, pooling = "partial", refresh = 0)
bg_mayo <- baggr(mayo, iter = 1e04)

plot(bg_mayo_partial)
baggr_compare("Mayo-Wilson" = bg_mayo, "Awasthi" = bg_a) %>% plot

# The difference in pooling metrics are not just a quirk of 
# baggr calculations, which take mean(SE)

rma.uni(yi = tau, sei = se, data = mayo)
rma.uni(yi = tau, sei = se, data = awasthi)
mean(pooling(bg_a)["mean",,1])

tau <- 0.37
mean(tau^2 / (tau^2 + bg_a$data$se^2))
tau^2 / (tau^2 + mean(bg_a$data$se^2))

tau <- 0.271
mean(tau^2 / (tau^2 + bg_a$data$se^2))
tau^2 / (tau^2 + mean(bg_a$data$se^2))

# inverse sampling variance method:
isq <- function(bg){
  se <- bg$data$se
  w <- 1/(se^2)
  k <- length(se)
  nu <- (k-1)*sum(w) / (sum(w)^2 - sum(w^2))
  tau <- treatment_effect(bg)[[2]]
  isq <- tau^2 / (tau^2 + nu)
  mint(isq)
}
isq(bg_a)
isq(bg_mayo)
