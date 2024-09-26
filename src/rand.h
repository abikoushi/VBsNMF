#include "RcppArmadillo.h"
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::mat rand_init(const arma::mat & alpha, const arma::rowvec & beta){
  arma::mat Z = alpha;
  for (int i=0; i<alpha.n_rows; i++) {
    for (int j=0; j<alpha.n_cols; j++) {
      Z(i,j) = arma::randg(arma::distr_param(alpha(i,j), 1/beta(j)));
    }
  }
  return Z;
}

/*
double NegativeSampling(const double & p1, 
                        const double & ns,
                        const arma::vec & v){
  const int n0 = R::rnbinom(ns, p1);
  const arma::vec prob = arma::ones<arma::vec>(v.n_rows);
  double out = 0;
  if(n0+ns < v.n_rows){
    out = sum(RcppArmadillo::sample(v, n0+ns, false, prob));    
  }else{
    out = sum(v);
  }
  //Rprintf("%f\n",out);
  return out;
}
*/
