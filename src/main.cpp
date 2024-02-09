#include "RcppArmadillo.h"
#define pi 3.141592653589793238462643383280

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
    
    if(iter == 100){
      Rcpp::stop("The result does not converge.");
    }
    
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
    
    if(iter == 100){
      Rcpp::stop("The result does not converge.");
    }
    
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
  double fisher = dat.n_rows/2.0;
  double x_new = xt + (dloglik(dat, xt) * (1/fisher));
  
  while(! conv_crit(xt, x_new, eps)){
    iter += 1;
    
    if(iter == 100){
      Rcpp::stop("The result does not converge.");
    }
    
    xt = x_new;
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
    
    if(iter == 100){
      Rcpp::stop("The result does not converge.");
    }
    
    x_prev = x_now;
    x_now = x_new;
    x_new = x_now - dloglik(dat, x_now) * ((x_now - x_prev)/(dloglik(dat, x_now) - dloglik(dat, x_prev)));
  }
  
  Rcpp::List result;
  result["xt"] = x_now;
  result["n_iter"] = iter;
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List optimQ3(arma::mat desMat, arma::vec Y, arma::vec b0, double eps){
  
  unsigned int iter = 1;
  arma::vec bt = b0;
  arma::vec xbeta = desMat * bt;
  arma::vec yhat = arma::exp(xbeta)/(1 + arma::exp(xbeta));;
  arma::mat W = arma::diagmat(arma::exp(xbeta)/arma::pow(1 + arma::exp(xbeta), 2.0));
  arma::mat invH = -arma::solve(desMat.t() * W * desMat, arma::eye(3, 3));
  arma::vec res = Y - yhat;
  arma::vec b_new = bt - (invH * desMat.t() * res);
  
  while(arma::norm(b_new - bt, 2) >= eps){
    iter += 1;
    
    if(iter == 100){
      Rcpp::stop("The result does not converge.");
    }
    
    bt = b_new;
    xbeta = desMat * bt;
    yhat = arma::exp(xbeta)/(1 + arma::exp(xbeta));;
    W = arma::diagmat(arma::exp(xbeta)/arma::pow(1 + arma::exp(xbeta), 2.0));
    invH = -arma::solve(desMat.t() * W * desMat, arma::eye(3, 3));
    res = Y - yhat;
    b_new = bt - (invH * desMat.t() * res);
  }

  Rcpp::List result;
  result["iter"] = iter;
  result["bt"] = bt;
  result["W"] = W;
  return result;
  
}
