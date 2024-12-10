#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

void readmtx(arma::uvec & row_i,
             arma::uvec & col_i,
             arma::vec & val,
             const std::string & readtxt,
             const arma::uvec & bag);

arma::umat randpick_c(int N1, int b_size);
