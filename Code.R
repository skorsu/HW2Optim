library(tidyverse)
library(Rcpp)
library(RcppArmadillo)

path <- "/Users/kevinkvp/Desktop/Github Repo/HW2Optim/src/"
path <- "/Users/kevin-imac/Desktop/Github - Repo/HW2Optim/src/"
sourceCpp(paste0(path, "main.cpp"))

### Q1
dat <- c(-8.86, -6.82, -4.03, -2.84, 0.14, 0.19, 0.24, 0.27, 0.49, 0.62, 0.76, 1.09,
         1.18, 1.32, 1.36, 1.58, 1.58, 1.78, 2.13, 2.15, 2.36, 4.05, 4.11, 4.12, 
         6.83)

sapply(seq(-10, 10, 0.01), ddloglik, x = dat) %>%
  plot(x = seq(-10, 10, 0.01), type = "l")

sapply(seq(-100, 100, 0.01), 
       function(x){dloglik(x = dat, theta = x)/ddloglik(x = dat, theta = x)}) %>%
  plot(x = seq(-100, 100, 0.01), type = "l")


10 - dloglik(x = dat, theta = 10)/ddloglik(x = dat, theta = 10)
33.87083 - dloglik(x = dat, theta = 33.87083)/ddloglik(x = dat, theta = 33.87083)


conv_crit(5, 4, 0.00001)
test <- bisect_q1(-10, 10, dat, eps = 1e-5)
nr_q1(x0 = 5, dat = dat, eps = 0.01)

dloglik(x = dat, theta = -10) * dloglik(x = dat, theta = 0)

uniroot(dloglik, interval = c(-10, 10), x = dat)

ddloglik(x = dat, theta = 10)

sum((-1 + (dat - 1)^2)/(1 + (dat - 1)^2)) * 2

1.5 - dloglik(x = dat, theta = 1.5)/ddloglik(x = dat, theta = 1.5)
