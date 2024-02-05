#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

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
  
  unsigned int iter = 0; // Iteration Counter 
  double xt = (a_init + b_init)/2.0; // Initialize x0
  double x_new = xt + eps + 1; // Initialize the x_new. Just a dummy one in order to let them pass the while criteria
  
  // Start the bisectioon method
  while(! conv_crit(xt, x_new, eps)){
    iter += 1;
    xt = (a_init + b_init)/2.0; // Propose new xt
    
    // Propose a new interval
    double eval_q = dloglik(dat, a_init);
    eval_q *= dloglik(dat, xt);
    
    if(eval_q <= 0){
      b_init = xt;
    } else {
      a_init = xt;
    }
  
    x_new = (a_init + b_init)/2.0; // Should we stop?
    
  }
  
  Rcpp::List result;
  result["a"] = a_init;
  result["b"] = b_init;
  result["xt"] = xt;
  result["n_iter"] = iter;
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List nr_q1(double x0, arma::vec dat, double eps){
  
  unsigned int iter = 0; // Iteration Counter
  
  // First iteration
  iter += 1;
  double xt = x0;
  double x_new = xt - (dloglik(dat, xt)/ddloglik(dat, xt));
  
  while(! conv_crit(xt, x_new, eps)){
    iter += 1;
    xt = x_new;
    x_new = xt - (dloglik(dat, xt)/ddloglik(dat, xt));
    std::cout << "xt: " << xt << std::endl;
    std::cout << "x_new: " << x_new << std::endl;
  }
  
  Rcpp::List result;
  result["xt"] = xt;
  result["iter"] = iter;
  return result;
  
}

