#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

void readmtx(arma::uvec & row_i,
             arma::uvec & col_i,
             arma::vec & val,
             const std::string & readtxt,
             const arma::uvec & bag);