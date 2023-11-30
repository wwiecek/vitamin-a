library(readr)
library(readxl)
library(tidyverse)

imdad2022 <- read_csv("input/imdad2022.csv", 
                      col_names = c("group", "tau", "se"))

mayo <- read_excel("input/mayo_wilson_fig3.xlsx") %>%
  mutate(tau = log(RR), 
         se = (log(uci) - log(lci)) / (2*1.96),
         se_alt = (log(uci) - log(RR)) / (1.96))
mayo$group[mayo$group == "Rahmatullah 1990"] <- "Rahmathullah 1990"
mayo$group[mayo$group == "Venkatarao 1996"] <- "Venkatrao 1996"
mayo$group[mayo$group == "Vilayaraghvan 1990"] <- "Vijayaraghavan 1990"


# Awasthi numbers are in Fig 4 here:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3647148/

# Note, that for Dibley 1996 there is something going on,
# probably this 0.01 is not precise enough
# ggplot(mayo, aes(x = group, y = 0, ymin = lci, ymax=uci)) + 
  # geom_point() + geom_errorbar() + coord_flip()

# I think LCI for Dibley study should be 0.0044
exp(-1.11 - 4.3)
# but maybe they didn't want to round to 0.00...

# You can see here how there is .01-.02 error in SE due to rounding off
# of LCI/UCI... but for Dibley it is .07



# Ross 1993 (survival) is also referred to as VAST 1993
awasthi <- data.frame(
  group = c("Sommer 1986", "Vijayaraghavan 1990", "Rahmathullah 1990", 
            "West 1991", "Daulaire 1992", "Herrera 1992", 
            "Arthur 1992", "Ross 1993 (survival)", "DEVTA 2013"),
  rr = c(0.66, 1.00, 0.46, 0.70, 0.74, 1.06, 0.30, 0.81, 0.96),
  se = c(0.202, 0.222, 0.220, 0.115, 0.150, 0.131, 0.466, 0.093, 0.036),
  lci = c(0.44, 0.65, 0.30, 0.56, 0.55, 0.82, 0.12, 0.68, 0.89),
  uci = c(0.97, 1.55, 0.71, 0.88, 0.99, 1.37, 0.75, 0.98, 1.03)
) %>% 
  # mutate(tau = log(rr)) %>% 
  mutate(tau = log(rr - se^2 / 2)) %>%
  mutate(lci2 = exp(tau - 1.96*se),
         uci2 = exp(tau + 1.96*se))
