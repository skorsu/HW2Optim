#include "RcppArmadillo.h"
#define pi 3.141592653589793238462643383280

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double loglik(arma::vec x, double theta){
  
  double result = 0.0;
  result -= arma::accu(arma::log(1 + arma::pow(x - theta, 2.0)));
  result -= (x.size() * std::log(pi));
  return result;
  
}

// [[Rcpp::export]]
double dloglik(arma::vec x, double theta){
  
  double result = 0.0;
  result += arma::accu((x - theta)/(1 + arma::pow(x - theta, 2)));
  result *= 2.0;
  return result;
  
}

// [[Rcpp::export]]
double ddloglik(arma::vec x, double theta){
  
  double result = 0.0;
  arma::vec xt2 = arma::pow(x - theta, 2.0);
  result += arma::accu((-1 + xt2)/arma::pow(1 + xt2, 2.0));
  result *= 2.0;
  return result;
  
}

// [[Rcpp::export]]
bool conv_crit(double x_old, double x_new, double eps){
  return (std::abs(x_new - x_old) < eps);
}

// [[Rcpp::export]]
Rcpp::List bisect_q1(double a_init, double b_init, arma::vec dat,
                     double eps){
  
  double a = a_init;
  double b = b_init;
  
  // First Iteration
  unsigned int iter = 1; // Iteration Counter 
  double xt = (a + b)/2.0;
  double eval_quan = dloglik(dat, a) * dloglik(dat, xt);
  if(eval_quan <= 0.0){
    b = xt;
  } else {
    a = xt;
  }
  double x_new = (a + b)/2.0;
  
  // Start the bisection method
  while(! conv_crit(xt, x_new, eps)){
    iter += 1;
    xt = x_new;
    eval_quan = dloglik(dat, a) * dloglik(dat, xt);
    if(eval_quan <= 0.0){
      b = xt;
    } else {
      a = xt;
    }
    x_new = (a + b)/2.0;
  }
  
  Rcpp::List result;
  result["xt"] = xt;
  result["n_iter"] = iter;
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List nr_q1(double x0, arma::vec dat, double eps){
  
  // First iteration
  unsigned int iter = 1; // Iteration Counter
  double xt = x0;
  double x_new = xt - (dloglik(dat, xt)/ddloglik(dat, xt));
  
  while(! conv_crit(xt, x_new, eps)){
    iter += 1;
    xt = x_new;
    x_new = xt - (dloglik(dat, xt)/ddloglik(dat, xt));
  }
  
  Rcpp::List result;
  result["xt"] = xt;
  result["n_iter"] = iter;
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List fs_q1(double x0, arma::vec dat, double eps){
  
  // First iteration
  unsigned int iter = 1;
  double xt = x0;
  double fisher = -ddloglik(dat, xt);
  double x_new = xt + (dloglik(dat, xt) * (1/fisher));
  
  while(! conv_crit(xt, x_new, eps)){
    iter += 1;
    xt = x_new;
    fisher = -ddloglik(dat, xt);
    x_new = xt + (dloglik(dat, xt) * (1/fisher));
  }
  
  Rcpp::List result;
  result["xt"] = xt;
  result["n_iter"] = iter;
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List sc_q1(double x0, double x1, arma::vec dat, double eps){
  
  // First iteration
  unsigned int iter = 1;
  double x_prev = x0;
  double x_now = x1;
  double x_new = x_now - (dloglik(dat, x_now) * ((x_now - x_prev)/(dloglik(dat, x_now) - dloglik(dat, x_prev))));
  
  while(! conv_crit(x_now, x_new, eps)){
    iter += 1;
    x_prev = x_now;
    x_now = x_new;
    x_new = x_now - dloglik(dat, x_now) * ((x_now - x_prev)/(dloglik(dat, x_now) - dloglik(dat, x_prev)));
  }
  
  Rcpp::List result;
  result["xt"] = x_now;
  result["n_iter"] = iter;
  return result;
  
}




































