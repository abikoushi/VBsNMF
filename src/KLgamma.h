#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat mat_digamma(arma::mat & a);

void up_log_gamma(arma::mat & logv, const arma::vec & a, const double & logb, const int & l);

double kld(const arma::mat & alpha_z,
           const arma::mat & beta_z,
           const arma::mat & alpha_w,
           const arma::mat & beta_w,
           const double & a,
           const double & b);
