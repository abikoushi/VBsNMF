#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "KLgamma.h"
#include "rand.h"

void up_log_gamma(arma::mat & logv, const arma::vec & a, const double & logb, const int & l){
  int K = logv.n_rows;
  for(int k=0;k<K;k++){
    logv(k,l) = R::digamma(a(k)) - logb;
  }
}

double kld(const arma::mat & alpha_z,
           const arma::mat & beta_z,
           const arma::mat & alpha_w,
           const arma::mat & beta_w,
           const double & a,
           const double & b){
  double lp = 0;
  for(int i=0; i<alpha_z.n_rows; i++){
    for(int l=0; l<alpha_z.n_cols; l++){
      lp += kl2gamma(a, b, alpha_z(i,l), beta_z(l));
    }
  }
  for(int i=0; i<alpha_w.n_rows; i++){
    for(int l=0; l<alpha_w.n_cols; l++){
      lp += kl2gamma(a, b, alpha_w(i,l), beta_w(l));
    }
  }
  return lp;
}

//shape parameters
double up_Az(arma::mat & alpha_z,
             const arma::mat & logZ,
             const arma::mat & logW,
             const arma::vec & y,
             const arma::uvec & rowi,
             const arma::uvec & coli,
             const double & a){
  //initialize by hyper parameter
  alpha_z.fill(a);
  double lp = 0;
  //inclement sufficient statistics
  for(int n=0; n<y.n_rows; n++){
    arma::rowvec r = exp(logZ.row(rowi(n)) + logW.row(coli(n)));
    double R = sum(r);
    r = y(n)*(r/R);
    alpha_z.row(rowi(n)) += r;
    lp +=  y(n)*log(R);
  }
  return lp;
}

double up_Aw(arma::mat & alpha_w,
             const arma::mat & logZ,
             const arma::mat & logW,
             const arma::vec & y,
             const arma::uvec & rowi,
             const arma::uvec & coli,
             const double & a){
  //initialize by hyper parameter
  alpha_w.fill(a);
  double lp = 0;
  //inclement sufficient statistics
  for(int n=0; n<y.n_rows; n++){
    arma::rowvec r = exp(logZ.row(rowi(n)) + logW.row(coli(n)));
    double R = sum(r);
    r = y(n)*(r/R);
    alpha_w.row(coli(n)) += r;
    lp +=  y(n)*log(R);
  }
  return lp;
}

/*
//rate parameters
double up_B(const arma::mat & alpha_z,
          const arma::mat & alpha_w,
          arma::mat & beta_z,
          arma::mat & beta_w,
          arma::mat & Z,
          arma::mat & W,
          arma::mat & logZ,
          arma::mat & logW,
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
    up_log_gamma(logW, alpha_w.col(l), log(beta_w(l)), l);
    //row z
    double B2 = sum(W.col(l));
    lp -= B2;
    beta_z.col(l) += B2;
    Z.col(l) = alpha_z.col(l)/beta_z(l);
    up_log_gamma(logZ, alpha_z.col(l), log(beta_z(l)), l);
  }
  return lp;
}
*/

double up_Bz(const arma::mat & alpha_z,
             arma::mat & beta_z,
             arma::mat & Z,
             arma::mat & W,
             arma::mat & logZ,
             const arma::uvec & rowi,
             const arma::uvec & coli,
             const double & b){
  int L = Z.n_cols;
  double lp = 0;
  beta_z.fill(b);
  //Rprintf("b\n");
  for(int l=0; l<L; l++){
    double B2 = sum(W.col(l));
    lp -= B2;
    beta_z.col(l) += B2;
    Z.col(l) = alpha_z.col(l)/beta_z(l);
    up_log_gamma(logZ, alpha_z.col(l), log(beta_z(l)), l);
  }
  return lp;
}

double up_Bw(const arma::mat & alpha_w,
             arma::mat & beta_w,
             arma::mat & Z,
             arma::mat & W,
             arma::mat & logW,
             const arma::uvec & rowi,
             const arma::uvec & coli,
             const double & b){
  int L = Z.n_cols;
  double lp = 0;
  beta_w.fill(b);
  for(int l=0; l<L; l++){
    double B1 = sum(Z.col(l)); 
    lp -= B1;
    beta_w.col(l) += B1;
    W.col(l) = alpha_w.col(l)/beta_w(l);
    up_log_gamma(logW, alpha_w.col(l), log(beta_w(l)), l);
    double B2 = sum(W.col(l));
    lp -= B2;
  }
  return lp;
}


double up_theta(arma::mat & alpha_z,
                arma::mat & alpha_w,
                arma::mat & beta_z,
                arma::mat & beta_w,
                arma::mat & Z,
                arma::mat & W,
                arma::mat & logZ,
                arma::mat & logW,
                const arma::vec & y,
                const arma::uvec & rowi,
                const arma::uvec & coli,
                const double & a,
                const double & b){
  int L = Z.n_cols;
  double lp = 0;
  up_Az(alpha_z,logZ,logW,y,rowi,coli,a);
  up_Bz(alpha_z,beta_z,Z,W,logZ,rowi,coli,b);
  lp += up_Aw(alpha_w,logZ,logW,y,rowi,coli,a);
  lp += up_Bw(alpha_w,beta_w,Z,W,logW,rowi,coli,b);
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
    double lp1 = up_theta(alpha_z, alpha_w, beta_z, beta_w, Z, W, logZ, logW, y, rowi, coli, a, b);
    lp(i) = lp1 + kld(alpha_z, beta_z, alpha_w, beta_w, a, b);
  }
  lp -= sum(lgamma(y+1)); 
  return List::create(Named("shape_row")=alpha_z,
                      Named("rate_row")=beta_z,
                      Named("shape_col")=alpha_w,
                      Named("rate_col")=beta_w,
                      Named("logprob")=lp);
}

////
//with NA
////

//shape parameters
double up_A(arma::mat & alpha_z,
            arma::mat & alpha_w,
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
    r = y(n)*(r/R);
    alpha_z.row(rowi(n)) += r;
    alpha_w.row(coli(n)) += r;
    lp +=  y(n)*log(R);
  }
  return lp;
}

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
    double lp_a = up_A(alpha_z, alpha_w, logZ, logW, y, rowi, coli, a);
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
/*
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
    r /= R;
    r *= y(n);
    alpha_z.row(rowi(n)) += r;
    alpha_w.row(coli(n)) += r;
  }
  return lp;
}
*/

/*
double NegativeSampling(const arma::vec & v){
  //future work sum(weight%v);
  double out = sum(v);
  return out;
}
*/

double up_Az_s(arma::mat & alpha_z,
               const arma::mat & logZ,
               const arma::mat & logW,
               const arma::vec & y,
               const arma::uvec & rowi,
               const arma::uvec & coli,
               const double & a){
  //initialize by hyper parameter
  //alpha_z.rows(uid_r).fill(a);
  alpha_z.fill(a);
  double lp = 0;
  //inclement sufficient statistics
  for(int n=0; n<y.n_rows; n++){
    int rn = rowi(n);
    int cn = coli(n);
    arma::rowvec r = exp(logZ.row(rn) + logW.row(cn));
    double R = sum(r);
    lp +=  y(n)*log(R) - lgamma(y(n)+1);
    r /= R;
    r *= y(n);
    alpha_z.row(rn) += r;
  }
  return lp;
}

double up_Aw_s(arma::mat & alpha_w,
               const arma::mat & logZ,
               const arma::mat & logW,
               const arma::vec & y,
               const arma::uvec & rowi,
               const arma::uvec & coli,
               const double & a){
  //initialize by hyper parameter
  //alpha_w.rows(uid_c).fill(a);
  alpha_w.fill(a);
  double lp = 0;
  //inclement sufficient statistics
  for(int n=0; n<y.n_rows; n++){
    int rn = rowi(n);
    int cn = coli(n);
    arma::rowvec r = exp(logZ.row(rn) + logW.row(cn));
    double R = sum(r);
    lp +=  y(n)*log(R) - lgamma(y(n)+1);
    r /= R;
    r *= y(n);
    alpha_w.row(cn) += r;
  }
  return lp;
}

double up_Bz_s(const arma::mat alpha_z,
               arma::mat & beta_z,
               arma::mat & Z,
               arma::mat & W,
               const arma::rowvec & SW,
               const double & weight, 
               const double & b){
  int L = Z.n_cols;
  double lp = 0;
  beta_z.fill(b);
  for(int l=0; l<L; l++){
    double B2=0;
    B2 += SW(l) +  weight*sum(W.col(l));
    lp -= B2;
    beta_z.col(l) += B2;
    Z.col(l) = alpha_z.col(l)/beta_z(l);
  }
  return lp;
}

double up_Bw_s(const arma::mat alpha_w,
               arma::mat & beta_w,
               arma::mat & Z,
               arma::mat & W,
               const arma::rowvec & SZ,
               const double & weight, 
               const double & b){
  int L = Z.n_cols;
  double lp = 0;
  beta_w.fill(b);
  for(int l=0; l<L; l++){
    double B1=0;
    B1 += SZ(l) + weight*(sum(Z.col(l)));
    lp -= B1;
    beta_w.col(l) += B1;
    W.col(l) = alpha_w.col(l)/beta_w(l);
  }
  return lp;
}

double up_theta_s(arma::mat & alpha_z,
                  arma::mat & alpha_w,
                  arma::mat & beta_z,
                  arma::mat & beta_w,
                  arma::mat & Z,
                  arma::mat & W,
                  arma::mat & logZ,
                  arma::mat & logW,
                  const arma::rowvec & SZ,
                  const arma::rowvec & SW,
                  const arma::vec & y,
                  const arma::uvec & rowi,
                  const arma::uvec & coli,
                  const double & a,
                  const double & b,
                  const double & weight){
  int L = Z.n_cols;
  double lp = 0;
  up_Az_s(alpha_z, logZ, logW, y, rowi, coli, a);
  up_Bz_s(alpha_z, beta_z, Z, W, SW, weight, b);
  logZ = mat_digamma(alpha_z).each_row() - log(beta_z);
  lp += up_Aw_s(alpha_w, logZ, logW, y, rowi, coli, a);
  lp += up_Bw_s(alpha_w, beta_w, Z, W, SZ, weight, b);
  logW = mat_digamma(alpha_w).each_row() - log(beta_w);
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
  arma::mat Z = randinit_gamma(alpha_z, beta_z);
  arma::mat W = randinit_gamma(alpha_w, beta_w);
  arma::mat logZ = log(Z);
  arma::mat logW = log(W);
  const double weight = (ns/N1);
  arma::rowvec SZ = beta_w - weight*sum(Z,0) - b;
  arma::rowvec SW = beta_z - weight*sum(W,0) - b;
  SZ.cols(find(SZ<0.0)).fill(0.0);
  SW.cols(find(SW<0.0)).fill(0.0);
  arma::vec lp = arma::zeros<arma::vec>(iter);
  for (int i=0; i<iter; i++) {
    double lp1 = up_theta_s(alpha_z, alpha_w, beta_z, beta_w, Z, W, logZ, logW, SZ, SW, y, rowi, coli, a, b, weight);
    lp(i) = lp1 + kld(alpha_z, beta_z, alpha_w, beta_w, a, b);
  }
  return List::create(Named("shape_row")=alpha_z,
                      Named("rate_row")=beta_z,
                      Named("shape_col")=alpha_w,
                      Named("rate_col")=beta_w,
                      Named("logprob")=lp);
}

////
//To Do
//doVB_pois_s_na : stochastic mini-batch with na
////
