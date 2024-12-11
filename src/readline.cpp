#include "readline.h"
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <string>
#include <vector>
using namespace Rcpp;

void readmtx(arma::uvec & row_i,
             arma::uvec & col_i,
             arma::vec & val,
             const std::string & readtxt,
             const arma::uvec & bag) {
  int x;
  int y;
  double v;
  std::ifstream file(readtxt);
  std::string str;    
  int index = 0;
  int n = 0;
  while (std::getline(file, str)){
      if(index == bag(n)+2){
        std::stringstream ss(str);
        std::vector<std::string> svec;
        while( ss.good() ){
          std::string substr;
          getline(ss, substr, ' ');
          svec.push_back(substr);
        }
        x = stoi(svec[0]);
        y = stoi(svec[1]);
        v = stod(svec[2]);
        row_i(n) = x-1;
        col_i(n) = y-1;
        val(n) = v;
        n++;
      }
      index++;
      if(n>=bag.n_rows){
        break;
      }
  }
}

// [[Rcpp::export]]
List read_mtx(const std::string & readtxt,
             const arma::uvec & bag){
  arma::uvec row_i(bag.n_rows);
  arma::uvec col_i(bag.n_rows);
  arma::vec val(bag.n_rows);
  readmtx(row_i,col_i,val, readtxt, bag);
  return List::create(row_i, col_i, val);
}

arma::umat randpick_c(int N1, int b_size){
  arma::uvec rind = arma::randperm(N1);
  int rem =  N1%b_size;
  int sb = N1/b_size;
  arma::umat subind;
  if(rem==0){
    subind = arma::zeros<arma::umat>(b_size, sb);
  }else{
    subind = arma::zeros<arma::umat>(b_size, sb+1);
  }
  if(rem==0){
    for(int i=0; i<sb; i++){
      subind.col(i) = rind.rows(b_size*i,b_size+b_size*i-1);
    }
  }else{
    for(int i=0; i<sb; i++){
      subind.col(i) = rind.rows(b_size*i, b_size+b_size*i-1);
    }
    arma::uvec r_rem  =  rind.rows(0, b_size-rem-1);
    subind.col(sb) = join_cols(rind.rows(b_size*sb, rind.n_rows-1), r_rem);
  }
  return subind;
}
