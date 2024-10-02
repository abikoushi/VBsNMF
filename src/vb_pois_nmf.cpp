#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "KLgamma.h"
#include "rand.h"

//shape parameters
double up_A(arma::mat & alpha_z,
            arma::mat & alpha_w,
            arma::mat & beta_z,
            arma::mat & beta_w,
            const arma::mat & logZ,
            const arma::mat & logW,
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
    arma::rowvec r = exp(logZ.row(rowi(n)) + logW.row(coli(n)));
    double R = sum(r);
    lp +=  y(n)*log(R) - lgamma(y(n)+1);
    r *= y(n)/R;
    alpha_z.row(rowi(n)) += r;
    alpha_w.row(coli(n)) += r;
  }
  return lp;
}

//rate parameters
double up_B(const arma::mat & alpha_z,
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
List doVB_pois(const arma::vec & y,
               const arma::uvec & rowi,
               const arma::uvec & coli,
               const int & Nr, const int & Nc, 
               const int & L,
               const int & iter,
               const double & a,
               const double & b){
  arma::mat Z = arma::randg<arma::mat>(Nr, L, arma::distr_param(a,1/b));
  arma::mat W = arma::randg<arma::mat>(Nc, L, arma::distr_param(a,1/b));
  arma::mat logZ = log(Z);
  arma::mat logW = log(W);
  arma::mat alpha_z = arma::ones<arma::mat>(Nr, L);
  arma::rowvec beta_z = arma::ones<arma::rowvec>(L);
  arma::mat alpha_w = arma::ones<arma::mat>(Nc, L);
  arma::rowvec beta_w = arma::ones<arma::rowvec>(L);
  arma::vec lp = arma::zeros<arma::vec>(iter);
  for (int i=0; i<iter; i++) {
    double lp_a = up_A(alpha_z, alpha_w, beta_z, beta_w, logZ, logW, y, rowi, coli, a);
    double lp_b = up_B(alpha_z, alpha_w, beta_z, beta_w, Z, W, rowi, coli, b);
    //lp.col(0).row(i) = lp_a ;
    //lp.col(1).row(i) = lp_b;
    //lp.col(2).row(i) = kld(alpha_z, beta_z, alpha_w, beta_w, a, b);
    lp(i) = lp_a+lp_b+kld(alpha_z, beta_z, alpha_w, beta_w, a, b);
    logZ = mat_digamma(alpha_z).each_row() - log(beta_z);
    logW = mat_digamma(alpha_w).each_row() - log(beta_w);
  }
  return List::create(Named("shape_row")=alpha_z,
                      Named("rate_row")=beta_z,
                      Named("shape_col")=alpha_w,
                      Named("rate_col")=beta_w,
                      Named("logprob")=lp);
}

////
//with NA
////
//rate parameters
double up_B_na(const arma::mat & alpha_z,
            const arma::mat & alpha_w,
            arma::mat & beta_z,
            arma::mat & beta_w,
            arma::mat & Z,
            arma::mat & W,
            const arma::uvec & rowi,
            const arma::uvec & coli,
            const arma::vec & wrow,
            const arma::vec & wcol,
            const double & b){
  int L = Z.n_cols;
  double lp = 0;
  beta_z.fill(b);
  beta_w.fill(b);
  //Rprintf("b\n");
  for(int l=0; l<L; l++){
    //col W
    double B1 = sum(wrow%Z.col(l)); 
    lp -= B1;
    beta_w.col(l) += B1;
    W.col(l) = alpha_w.col(l)/beta_w(l);
    //row z
    double B2 =sum(wcol%W.col(l));
    lp -= B2;
    beta_z.col(l) += B2;
    Z.col(l) = alpha_z.col(l)/beta_z(l);
  }
  return lp;
}

// [[Rcpp::export]]
List doVB_pois_na(const arma::vec & y,
               const arma::uvec & rowi,
               const arma::uvec & coli,
               const int & Nr, const int & Nc,
               const int & L,
               const int & iter,
               const double & a,
               const double & b,
               const arma::vec & wrow,
               const arma::vec & wcol){
  arma::mat Z = arma::randg<arma::mat>(Nr, L);
  arma::mat W = arma::randg<arma::mat>(Nc, L);
  arma::mat logZ = log(Z);
  arma::mat logW = log(W);
  arma::mat alpha_z = arma::ones<arma::mat>(Nr, L);
  arma::rowvec beta_z = arma::ones<arma::rowvec>(L);
  arma::mat alpha_w = arma::ones<arma::mat>(Nc, L);
  arma::rowvec beta_w = arma::ones<arma::rowvec>(L);
  arma::vec lp = arma::zeros<arma::vec>(iter);
  for (int i=0; i<iter; i++) {
    double lp_a = up_A(alpha_z, alpha_w, beta_z, beta_w, logZ, logW, y, rowi, coli, a);
    double lp_b = up_B_na(alpha_z, alpha_w, beta_z, beta_w, Z, W, rowi, coli, wrow, wcol, b);
    lp(i) = lp_a+lp_b+kld(alpha_z, beta_z, alpha_w, beta_w, a, b);
    logZ = mat_digamma(alpha_z).each_row() - log(beta_z);
    logW = mat_digamma(alpha_w).each_row() - log(beta_w);
  }
  return List::create(Named("shape_row")=alpha_z,
                      Named("rate_row")=beta_z,
                      Named("shape_col")=alpha_w,
                      Named("rate_col")=beta_w,
                      Named("logprob")=lp);
}

////
//s : stochastic mini-batch
////
//shape parameters
double up_A_s(arma::mat & alpha_z,
            arma::mat & alpha_w,
            arma::mat & beta_z,
            arma::mat & beta_w,
            const arma::mat & logZ,
            const arma::mat & logW,
            const arma::vec & y,
            const arma::uvec & rowi,
            const arma::uvec & coli,
            const arma::uvec & uid_r,
            const arma::uvec & uid_c,
            const double & a){
  //initialize by hyper parameter
  alpha_z.rows(uid_r).fill(a);
  alpha_w.rows(uid_c).fill(a);
  double lp = 0;
  //inclement sufficient statistics
  for(int n=0; n<y.n_rows; n++){
    arma::rowvec r = exp(logZ.row(rowi(n)) + logW.row(coli(n)));
    double R = sum(r);
    lp +=  y(n)*log(R) - lgamma(y(n)+1);
    r *= y(n)/R;
    alpha_z.row(rowi(n)) += r;
    alpha_w.row(coli(n)) += r;
  }
  return lp;
}

double NegativeSampling(const arma::vec & v){
  //future work sum(weight%v);
  double out = sum(v);
  return out;
}

//rate parameters
double up_B_s(const arma::mat alpha_z,
              const arma::mat alpha_w,
              arma::mat & beta_z,
              arma::mat & beta_w,
              arma::mat & Z,
              arma::mat & W,
              const double & weight, 
              const arma::uvec & rowi,
              const arma::uvec & coli,
              const double & b){
  int L = Z.n_cols;
  double lp = 0;
  beta_z.fill(b);
  beta_w.fill(b);
  for(int l=0; l<L; l++){
    //col W
    double B1=0;
    // for(int i=0; i < rowi.n_rows; i++){
    //   B1 += Z(rowi(i),l);
    // }
    B1 += weight*NegativeSampling(Z.col(l));
    lp -= B1;
    beta_w.col(l) += B1;
    W.col(l) = alpha_w.col(l)/beta_w(l);
    //row z
    double B2=0;
    // for(int i=0; i < coli.n_rows; i++){
    //   B2 += W(coli(i),l);
    // }
    B2 += weight*NegativeSampling(W.col(l));
    lp -= B2;
    beta_z.col(l) += B2;
    Z.col(l) = alpha_z.col(l)/beta_z(l);
  }
  return lp;
}


// [[Rcpp::export]]
List doVB_pois_s(const arma::vec & y,
                 const arma::uvec & rowi,
                 const arma::uvec & coli,
                 const int & L,
                 const int & iter,
                 const double & a,
                 const double & b,
                 const double & N1, 
                 arma::mat alpha_z, 
                 arma::rowvec beta_z,
                 arma::mat alpha_w, 
                 arma::rowvec beta_w){
  const double ns = y.n_rows;
  const arma::uvec uid_r = unique(rowi); //for a
  const arma::uvec uid_c = unique(coli); //for a
  arma::mat Z = rand_init(alpha_z, beta_z);
  arma::mat W = rand_init(alpha_w, beta_w);
  arma::mat logZ = log(Z);
  arma::mat logW = log(W);
  const double weight = (ns/N1);
  arma::vec lp = arma::zeros<arma::vec>(iter);
  for (int i=0; i<iter; i++) {
    double lp_a = up_A_s(alpha_z, alpha_w, beta_z, beta_w, logZ, logW, y, rowi, coli, uid_r, uid_c, a);
    //Rprintf("a\n");
    double lp_b = up_B_s(alpha_z, alpha_w, beta_z, beta_w, Z, W, weight, rowi, coli, b);
    //Rprintf("b\n");
    lp(i) = lp_a+lp_b+kld(alpha_z, beta_z, alpha_w, beta_w, a, b);
    logZ = mat_digamma(alpha_z).each_row() - log(beta_z);
    logW = mat_digamma(alpha_w).each_row() - log(beta_w);
  }
  return List::create(Named("shape_row")=alpha_z,
                      Named("rate_row")=beta_z,
                      Named("shape_col")=alpha_w,
                      Named("rate_col")=beta_w,
                      Named("logprob")=lp);
}