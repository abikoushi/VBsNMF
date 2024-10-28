#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "KLgamma.h"
#include "rand.h"


//shape parameters
double up_Aem(arma::mat & alpha_z,
            arma::mat & alpha_w,
            const arma::mat & Z,
            const arma::mat & W,
            const arma::vec & y,
            const arma::uvec & rowi,
            const arma::uvec & coli,
            const double & a){
  //initialize by hyper parameter
  alpha_z.fill(a);
  alpha_w.fill(a);
  double lp = 0;
  //inclement sufficient statistics
  for(int n=0; n<y.n_rows; n++){
    arma::rowvec r = Z.row(rowi(n)) % W.row(coli(n));
    double R = sum(r);
    r = y(n)*(r/R);
    alpha_z.row(rowi(n)) += r;
    alpha_w.row(coli(n)) += r;
    lp +=  y(n)*log(R);
  }
  return lp;
}

//rate parameters
double up_Bem(const arma::mat & alpha_z,
            const arma::mat & alpha_w,
            arma::mat & beta_z,
            arma::mat & beta_w,
            arma::mat & Z,
            arma::mat & W,
            const arma::uvec & rowi,
            const arma::uvec & coli,
            const double & b){
  int L = Z.n_cols;
  double lp = 0;
  beta_z.fill(b);
  beta_w.fill(b);
  //Rprintf("b\n");
  for(int l=0; l<L; l++){
    //col W
    double B1 = sum(Z.col(l)); 
    lp -= B1;
    beta_w.col(l) += B1;
    W.col(l) = alpha_w.col(l)/beta_w(l);
    //row z
    double B2 = sum(W.col(l));
    lp -= B2;
    beta_z.col(l) += B2;
    Z.col(l) = alpha_z.col(l)/beta_z(l);
  }
  return lp;
}


// [[Rcpp::export]]
List doEM_pois(const arma::vec & y,
               const arma::uvec & rowi,
               const arma::uvec & coli,
               const int & Nr, const int & Nc, 
               const int & L,
               const int & iter,
               const double & a,
               const double & b){
  arma::mat Z = arma::randg<arma::mat>(Nr, L);
  arma::mat W = arma::randg<arma::mat>(Nc, L);
  arma::vec lp = arma::zeros<arma::vec>(iter);
  for (int i=0; i<iter; i++) {
    arma::mat alpha_z = arma::zeros<arma::mat>(Nr, L);
    arma::rowvec beta_z = arma::zeros<arma::rowvec>(L);
    arma::mat alpha_w = arma::zeros<arma::mat>(Nc, L);
    arma::rowvec beta_w = arma::zeros<arma::rowvec>(L);
    double lp_a = up_Aem(alpha_z, alpha_w, Z, W, y, rowi, coli, a);
    double lp_b = up_Bem(alpha_z, alpha_w, beta_z, beta_w, Z, W, rowi, coli, b);
    lp(i) = lp_a+lp_b;
  }
  lp -= sum(lgamma(y+1)); 
  return List::create(Named("Z")=Z,
                      Named("W")=W,
                      Named("logprob")=lp);
}
