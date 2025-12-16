#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "KLgamma.h"
#include "rand.h"
#include "readline.h"
#include "lr.h"
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppProgress)]]

// [[Rcpp::export]]
double lr_default(const double & t,
                  const double & delay,
                  const double & forgetting){
  return pow(t+delay, -forgetting);
}

//shape parameters
double up_A(arma::mat & alpha,
             const arma::mat & logV_up,
             const arma::mat & logV_con,
             const arma::vec & y,
             const arma::uvec & up_i,
             const arma::uvec & con_i,
             const double & a){
  //initialize by hyper parameter
  alpha.fill(a);
  double lp = 0;
  //inclement sufficient statistics
  arma::mat r = exp(logV_up.rows(up_i) + logV_con.rows(con_i));
  arma::vec R = sum(r, 1);
  r.each_col() /= R;
  r.each_col() %= y;
  alpha.rows(up_i) += r;
  lp +=  sum( y % log(R) );
  return lp;
}

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

double up_B(const arma::mat & alpha,
             arma::mat & beta,
             arma::mat & V_up,
             arma::mat & V_con,
             arma::mat & logV_up,
             const double & b){
  int L = V_up.n_cols;
  double lp = 0;
  beta.fill(b);
  //Rprintf("b\n");
  for(int l=0; l<L; l++){
    double B2 = sum(V_con.col(l));
    lp -= B2;
    beta.col(l) += B2;
    V_up.col(l) = alpha.col(l)/beta(l);
    up_log_gamma(logV_up, alpha.col(l), log(beta(l)), l);
  }
  return lp;
}

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
  double lp = 0;
  up_A(alpha_z, logZ, logW, y, rowi, coli, a);
  up_B(alpha_z, beta_z, Z, W, logZ, b);
  lp += up_A(alpha_w, logW, logZ, y, coli, rowi, a);
  lp += up_B(alpha_w, beta_w, W, Z, logW, b);
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

double up_Az_s(arma::mat & alpha_z,
               const arma::mat & logZ,
               const arma::mat & logW,
               const arma::vec & y,
               const arma::uvec & rowi,
               const arma::uvec & coli,
               const double & a){
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
    //alpha_z.row(rn).print();
    alpha_z.row(rn) += r;
    //Rprintf("n\n");
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
  double lp = 0;

  up_Az_s(alpha_z, logZ, logW, y, rowi, coli, a);
  up_Bz_s(alpha_z, beta_z, Z, W, SW, weight, b);
  logZ = mat_digamma(alpha_z).each_row() - log(beta_z);
  lp += up_Aw_s(alpha_w, logZ, logW, y, rowi, coli, a);
  lp += up_Bw_s(alpha_w, beta_w, Z, W, SZ, weight, b);
  logW = mat_digamma(alpha_w).each_row() - log(beta_w);
  return lp;
}

double doVB_pois_s_sub(const arma::vec & y,
                 const arma::uvec & rowi,
                 const arma::uvec & coli,
                 const int & L,
                 const int & iter,
                 const double & a,
                 const double & b,
                 const double & N1, 
                 arma::mat & alpha_z, 
                 arma::rowvec & beta_z,
                 arma::mat & alpha_w, 
                 arma::rowvec & beta_w){
  arma::mat Z = randinit_gamma(alpha_z, beta_z);
  arma::mat W = randinit_gamma(alpha_w, beta_w);
  arma::mat logZ = log(Z);
  arma::mat logW = log(W);
  const double ns = y.n_rows;
  const double weight = (ns/N1);
  arma::rowvec SZ = beta_w - weight*sum(Z, 0) - b;
  arma::rowvec SW = beta_z - weight*sum(W, 0) - b;
  SZ.cols(find(SZ<0.0)).fill(0.0);
  SW.cols(find(SW<0.0)).fill(0.0);
  double lp = 0.0;
  for (int i=0; i<iter; i++) {
    double lp1 = up_theta_s(alpha_z, alpha_w, beta_z, beta_w, Z, W, logZ, logW, SZ, SW, y, rowi, coli, a, b, weight);
    lp = lp1 + kld(alpha_z, beta_z, alpha_w, beta_w, a, b);
  }
  return lp;
}


// [[Rcpp::export]]
List doVB_pois_s_mtx(const std::string & file_path,
                 const int & L,
                 const int & iter,
                 const int & subiter,
                 const double & a,
                 const double & b,
                 const double & N1,
                 const int & Nr, const int & Nc,
                 const int & ns,
                 const arma::vec & lr_param,
                 const std::string & lr_type,
                 const bool & display_progress){
  arma::mat Z = arma::randg<arma::mat>(Nr, L);
  arma::mat W = arma::randg<arma::mat>(Nc, L);
  const double weight = (ns/N1);
  arma::mat alpha_z = arma::ones<arma::mat>(Nr, L);
  arma::rowvec beta_z = arma::ones<arma::rowvec>(L);
  arma::mat alpha_w = arma::ones<arma::mat>(Nc, L);
  arma::rowvec beta_w = arma::ones<arma::rowvec>(L);
  arma::vec lp = arma::zeros<arma::vec>(iter);
  std::unique_ptr<lr> g;
  set_lr_method(g, lr_type);
  Progress pb(iter, display_progress);
  for(int epoc=0; epoc<iter; epoc++){
    arma::uvec row_i(ns);
    arma::uvec col_i(ns);
    arma::vec val(ns);
    arma::umat bags = randpick_c(N1, ns);
    for(int step = 0; step < bags.n_cols; step++){
      arma::uvec bag = sort(bags.col(step));
      readmtx(row_i, col_i, val, file_path, bag);
      
      arma::uvec uid_r = unique(row_i);
      arma::uvec uid_c = unique(col_i);
      
      arma::mat alpha_zs = alpha_z.rows(uid_r);
      arma::mat alpha_ws = alpha_w.rows(uid_c);
      arma::rowvec beta_zs = beta_z;
      arma::rowvec beta_ws = beta_w;
      
      rankindex(row_i, uid_r);
      rankindex(col_i, uid_c);

      lp(epoc) += doVB_pois_s_sub(val, row_i, col_i, L, subiter, a, b, 
         N1, alpha_zs, beta_zs, alpha_ws, beta_ws);
      
      double rho = g -> lr_t(epoc, lr_param);
      double rho2 = 1-rho;
      alpha_z.rows(uid_r) = rho2*alpha_z.rows(uid_r) + rho*alpha_zs;
      beta_z = rho2*beta_z + rho*beta_zs;
      alpha_w.rows(uid_c) = rho2*alpha_w.rows(uid_c)+ rho*alpha_ws;
      beta_w = rho2*beta_w + rho*beta_ws;
    }
    pb.increment();
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
