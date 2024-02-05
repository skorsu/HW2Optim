#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double dloglik(arma::vec x, double theta){
  
  double result = 0.0;
  result += arma::accu((x - theta)/(1 + arma::pow(x - theta, 2)));
  result *= 2.0;
  return result;
  
}