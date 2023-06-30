library(readxl)
library(tidyverse)
library(baggr)
library(metafor)

set.seed(1990)
theme_set(theme_minimal(base_size = 10))

mayo <- read_excel("input/mayo_wilson_fig3.xlsx") %>%
  mutate(tau = log(RR), 
         se = (log(uci) - log(lci)) / (2*1.96),
         se_alt = (log(uci) - log(RR)) / (1.96))

awasthi <- data.frame(
  group = c("Sommer 1986", "Vijayaraghavan 1990", "Rahmathullah 1990", 
            "West 1990", "Daulaire 1992", "Herrera 1992", 
            "Arthur 1992", "VAST 1993", "DEVTA 2013"),
  rr = c(0.66, 1.00, 0.46, 0.70, 0.74, 1.06, 0.30, 0.81, 0.96),
  se = c(0.202, 0.222, 0.220, 0.115, 0.150, 0.131, 0.466, 0.093, 0.036),
  lci = c(0.44, 0.65, 0.30, 0.56, 0.55, 0.82, 0.12, 0.68, 0.89),
  uci = c(0.97, 1.55, 0.71, 0.88, 0.99, 1.37, 0.75, 0.98, 1.03)
) %>% 
  # mutate(tau = log(rr)) %>% 
  mutate(tau = log(rr - se^2 / 2)) %>%
  mutate(lci2 = exp(tau - 1.96*se),
         uci2 = exp(tau + 1.96*se))

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
