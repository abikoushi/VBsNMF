#include "RcppArmadillo.h"
#include "KLgamma.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double KLnorm(double mq, double bq, double tau){
  double out = 0.5*(tau/bq - 1.0 - log(tau) + log(bq));
  return isnan(out)?0:out;
}

double KLnorm_mat(arma::mat Z, arma::mat W,
                  arma::mat prec_z, arma::mat prec_w,
                  double prior_prec){
  double klv = 0;
  for(int l=0; l<Z.n_cols; l++){
    for(int i=0; i<Z.n_rows; i++){
      klv += KLnorm(Z(i,l), prec_z(i,l), prior_prec);
    }
    for(int i=0; i<W.n_rows; i++){
      klv += KLnorm(W(i,l), prec_w(i,l), prior_prec);
    }
  }
  return klv;
}

//arrow negative value
double up_theta(const arma::vec & y,
                const arma::uvec & rowi,
                const arma::uvec & coli,
                  arma::mat & Z,
                  arma::mat & W,
                  arma::mat & Z2,
                  arma::mat & W2,
                  arma::mat & prec_z,
                  arma::mat & prec_w,
                  const double & obs_prec,
                  const double & prior_prec) {
  double lp = 0;
  for(int l=0; l<Z.n_cols; l++){
    arma::vec num_z = Z.col(l);
    num_z.fill(0);
    arma::vec num_w = W.col(l);
    num_w.fill(0);
    arma::vec den_z = Z.col(l);
    den_z.fill(0);
    arma::vec den_w = W.col(l);
    den_w.fill(0);
  for(int n=0; n<y.n_rows; n++){
      double zw = arma::dot(Z.row(rowi(n)), W.row(coli(n)));
    lp += pow(y(n)-zw,2);
    double resid = y(n) - zw/(Z(rowi(n),l)*W(coli(n),l));
        num_z(rowi(n)) += W(coli(n),l)*resid;
        num_w(coli(n)) += Z(rowi(n),l)*resid;
        den_z(rowi(n)) += W2(coli(n),l);
        den_w(coli(n)) += Z2(rowi(n),l);
    }
  Z.col(l) = num_z/(den_z + obs_prec/prior_prec);
  prec_z.col(l) = obs_prec*den_z + prior_prec;
  W.col(l) = num_w/(den_w + obs_prec/prior_prec);
  prec_w.col(l) = obs_prec*den_w + prior_prec;
  Z2.col(l) = pow(Z.col(l),2) + 1.0/prec_z.col(l);
  W2.col(l) = pow(W.col(l),2) + 1.0/prec_w.col(l);
  }
  return lp;
}

double up_lambda(double & lambda,
                 const double & ahat, 
                 const double & a, const double & b){
  double tmp = 0;
  double bhat = tmp+b;
  lambda = ahat/bhat;
  return (-lambda*tmp) + 0.5*std::log(lambda) - kl2gamma(a, b, ahat, bhat);
}

// [[Rcpp::export]]
List doVB_norm(const arma::vec & y,
               const arma::uvec & rowi,
               const arma::uvec & coli,
               const int & Nr, const int & Nc, 
               const int & L,
               const int & iter,
               const double & prior_prec,
               const double & a,
               const double & b){
  arma::mat Z = arma::randn<arma::mat>(Nr, L);
  arma::mat W = arma::randn<arma::mat>(Nc, L);
  arma::mat Z2 = pow(Z, 2);
  arma::mat W2 = pow(W, 2);
  arma::mat prec_z = Z2;
  arma::mat prec_w = W2;
  double obs_prec = 1.0;
  //double obs_prec = a/b;
  arma::vec lp = arma::zeros<arma::vec>(iter);
  double N = Nr*Nc;
  double ahat = 0.5*N+a;
  for (int i=0; i<iter; i++) {
    lp(i) = up_theta(y, rowi, coli, Z, W, Z2, W2, prec_z, prec_w, obs_prec, prior_prec);
    //up_lambda(obs_prec, ahat, a, b);
  }
  return List::create(Named("mean_z") = Z, Named("mean_w") = W,
                      Named("prec_z") = prec_z, Named("prec_w") = prec_w,
                      Named("lp")=lp);
}
