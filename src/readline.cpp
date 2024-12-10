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
  for(int n = 0 ; n < bag.n_elem; n++){
    int x;
    int y;
    double v;
    //readaline_mtx(x, y, v, readtxt, bag(n)+2);
    
    std::ifstream file(readtxt);
    std::string str;
    
    int index = 0;
    while (std::getline(file, str))
    {
      if(index == bag(n)+2){
        std::stringstream ss(str);
        std::vector<std::string> svec;
        while( ss.good() ){
          std::string substr;
          getline(ss, substr, ' ');
          svec.push_back(substr);
        }
        int K = svec.size();
        x = stoi(svec[0]);
        y = stoi(svec[1]);
        v = stod(svec[2]);
        row_i(n) = x-1;
        col_i(n) = y-1;
        val(n) = v;
      }
      index++;
      if(index>bag(bag.n_rows-1)){
        break;
      }
    }
  }
  //return List::create(row_i, col_i, val);
}
