#include "RcppArmadillo.h"
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double klgamma_sub(double a, double b, double c, double d){
  return - c / (a*d) - log(a) / b + lgamma(b) + (R::digamma(d) + log(c))/(b-1);
}

double kl2gamma(double a1, double b1, double a2, double b2){
  return klgamma_sub(a2,b2,a2,b2) - klgamma_sub(a1,b1,a2,b2);
}

arma::mat mat_digamma(arma::mat a){
  int K = a.n_rows;
  int L = a.n_cols;
  arma::mat out(K,L);
  for(int k=0;k<K;k++){
    for(int l=0;l<L;l++){
      out(k,l) = R::digamma(a(k,l));
    }
  }
  return out;
}

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

double lowerbound_logML_pois(const arma::vec & y,
                             const arma::mat & alpha_z,
                             const arma::mat & beta_z,
                             const arma::mat & alpha_w,
                             const arma::mat & beta_w,
                             const arma::mat & Z,
                             const arma::mat & W,
                             const arma::mat & logZ,
                             const arma::mat & logW,
                             const double & a,
                             const double & b){
  double lp = 0;
  lp += sum(sum((a-1)*logZ - b*Z,0) +a*log(beta_z) - Z.n_rows*std::lgamma(a));
  lp -= accu((alpha_z-1)%logZ - Z.each_row()%beta_z + alpha_z.each_row()%log(beta_z) - lgamma(alpha_z));
  lp += sum(sum((a-1)*logW - b*W,0) + a*log(beta_w) - W.n_rows*W.n_elem*std::lgamma(a));
  lp -= accu((alpha_w-1)%logW - W.each_row()%beta_w + alpha_w.each_row()%log(beta_w) - lgamma(alpha_w));
  return lp;
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

arma::mat rand_init(const arma::mat & alpha, const arma::rowvec & beta){
  arma::mat Z = alpha;
  for (int i=0; i<alpha.n_rows; i++) {
    for (int j=0; j<alpha.n_cols; j++) {
      Z(i,j) = arma::randg(arma::distr_param(alpha(i,j), 1/beta(j)));
    }
  }
  return Z;
}

// [[Rcpp::export]]
List doVB_pois(const arma::vec & y,
               const arma::uvec & rowi,
               const arma::uvec & coli,
               const int & N, const int & M,
               const int & L,
               const int & iter,
               const double & a,
               const double & b,
               arma::mat & alpha_z, arma::rowvec & beta_z,
               arma::mat & alpha_w, arma::rowvec & beta_w){
  arma::mat Z = rand_init(alpha_z, beta_z);
  arma::mat W = rand_init(alpha_w, beta_w);
  arma::mat logZ = log(Z);
  arma::mat logW = log(W);
  arma::mat lp(iter,3);
  for (int i=0; i<iter; i++) {
    double lp_a = up_A(alpha_z, alpha_w, beta_z, beta_w, logZ, logW, y, rowi, coli, a);
    double lp_b = up_B(alpha_z, alpha_w, beta_z, beta_w, Z, W, rowi, coli, b);
    lp.col(0).row(i) = lp_a ;
    lp.col(1).row(i) = lp_b;
    lp.col(2).row(i) = kld(alpha_z, beta_z, alpha_w, beta_w, a, b);
    logZ = mat_digamma(alpha_z).each_row() - log(beta_z);
    logW = mat_digamma(alpha_w).each_row() - log(beta_w);
  }
  return List::create(Named("shape_row")=alpha_z,
                      Named("rate_row")=beta_z,
                      Named("shape_col")=alpha_w,
                      Named("rate_col")=beta_w,
                      Named("logprob")=lp);
}