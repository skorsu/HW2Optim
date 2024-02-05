library(tidyverse)
library(Rcpp)
library(RcppArmadillo)

path <- "/Users/kevin-imac/Desktop/Github - Repo/HW2Optim/src/"
sourceCpp(paste0(path, "main.cpp"))

### Q1
dat <- c(-8.86, -6.82, -4.03, -2.84, 0.14, 0.19, 0.24, 0.27, 0.49, 0.62, 0.76, 1.09,
         1.18, 1.32, 1.36, 1.58, 1.58, 1.78, 2.13, 2.15, 2.36, 4.05, 4.11, 4.12, 
         6.83)

sapply(seq(-10, 10, 0.01), dloglik, x = dat) %>%
  plot(type = "l")
