library(tidyverse)
library(Rcpp)
library(RcppArmadillo)

# path <- "/Users/kevinkvp/Desktop/Github Repo/HW2Optim/src/"
path <- "/Users/kevin-imac/Desktop/Github - Repo/HW2Optim/src/"
sourceCpp(paste0(path, "main.cpp"))

### Q1
dat <- c(-8.86, -6.82, -4.03, -2.84, 0.14, 0.19, 0.24, 0.27, 0.49, 0.62, 0.76, 1.09,
         1.18, 1.32, 1.36, 1.58, 1.58, 1.78, 2.13, 2.15, 2.36, 4.05, 4.11, 4.12, 
         6.83)
add_dat <- c(-8.34, -1.73, -0.40, -0.24, 0.60, 0.94, 1.05, 1.06, 1.45, 1.50, 
             1.54, 1.72, 1.74, 1.88, 2.04, 2.16, 2.39, 3.01, 3.01, 3.08, 4.66,
             4.99, 6.01, 7.06, 25.45)

sapply(seq(-100, 100, 0.01), loglik, x = c(dat)) %>%
  plot(x = seq(-100, 100, 0.01), type = "l")


median(c(dat, add_dat))

sapply(seq(-10, 10, 0.01), ddloglik, x = c(dat, add_dat)) %>%
  plot(x = seq(-10, 10, 0.01), type = "l")
abline(h = 0)
abline(v = 0)
abline(v = 2.5)

sapply(seq(0, 10, 0.01), 
       function(x){dloglik(x = dat, theta = x)}) %>%
  plot(x = seq(0, 10, 0.01), type = "l")


10 - dloglik(x = dat, theta = 10)/ddloglik(x = dat, theta = 10)
33.87083 - dloglik(x = dat, theta = 33.87083)/ddloglik(x = dat, theta = 33.87083)


conv_crit(5, 4, 0.00001)
bisect_q1(-10000, 10000, dat, eps = 1e-5)
nr_q1(x0 = 0.1, dat = c(dat, add_dat), eps = 1e-5)
fs_q1(x0 = 1e-10, dat = dat, eps = 1e-5)
sc_q1(x0 = 0, x1 = 0.5, dat = c(dat, add_dat), eps = 1e-5)

dloglik(x = dat, theta = -10) * dloglik(x = dat, theta = 0)

uniroot(dloglik, interval = c(-10, 10), x = dat)

uniroot(ddloglik, interval = c(0, 2.5), x = dat)

ddloglik(x = dat, theta = 10)

sum((-1 + (dat - 1)^2)/(1 + (dat - 1)^2)) * 2

1.5 - dloglik(x = dat, theta = 1.5)/ddloglik(x = dat, theta = 1.5)

apply(data.frame(expand.grid(seq(-1, 1, 0.01), seq(-1, 1, 0.01))), 1, 
      function(x){dloglik(dat, x[1]) - dloglik(dat, x[2])}) %>%
  hist()


sapply(seq(-10, 10, 0.01), ddloglik, x = c(dat, add_dat)) %>%
  plot(x = seq(-10, 10, 0.01), type = "l")
