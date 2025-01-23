#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace Rcpp;

///
//read
///
// [[Rcpp::export]]
arma::mat readRowsFromBinary(std::string filepath, arma::uvec rows, int ncols) {
  std::ifstream in(filepath, std::ios::binary);
  if (!in.is_open()) {
    Rcpp::stop("Failed to open file for reading.");
  }
  int nrows = rows.n_elem;
  arma::mat result(nrows, ncols);
  
  for (arma::uword i = 0; i < nrows; ++i) {
    size_t offset = rows(i) * ncols * sizeof(double);
    in.seekg(offset, std::ios::beg);
    
    if (!in) {
      Rcpp::stop("Error seeking to the specified position in the file.");
    }
    
    std::vector<double> buffer(ncols);
    in.read(reinterpret_cast<char*>(buffer.data()), ncols * sizeof(double));
    
    if (!in) {
      Rcpp::stop("Error reading from the file.");
    }
    result.row(i) = arma::rowvec(buffer);
  }
  
  in.close();
  return result;
}

// [[Rcpp::export]]
arma::umat readRowsFromBinary_umat(std::string filepath, arma::uvec rows, int ncols) {
  std::ifstream in(filepath, std::ios::binary);
  if (!in.is_open()) {
    Rcpp::stop("Failed to open file for reading.");
  }
  
  // 読み込む行数を取得し、結果格納用の行列を初期化
  int nrows = rows.n_elem;
  arma::umat result(nrows, ncols);
  
  // 各行を読み込む
  for (arma::uword i = 0; i < nrows; ++i) {
    size_t offset = rows[i] * ncols * sizeof(unsigned int);
    in.seekg(offset, std::ios::beg);
    
    if (!in) {
      Rcpp::stop("Error seeking to the specified position in the file.");
    }
    
    std::vector<unsigned int> buffer(ncols);
    in.read(reinterpret_cast<char*>(buffer.data()), ncols * sizeof(unsigned int));
    
    if (!in) {
      Rcpp::stop("Error reading from the file.");
    }
    
    // 読み込んだデータを結果の行列に格納
    for (arma::uword j = 0; j < ncols; ++j) {
      result(i, j) = buffer[j];
    }
  }
  
  in.close();
  return result;
}

// [[Rcpp::export]]
List read_bin(const std::string & filepath_x,
              const std::string & filepath_y,
              const arma::uvec & bag,
              const int & n_dim){
  arma::umat X(bag.n_rows, n_dim-1);
  arma::vec val(bag.n_rows);
  //readBin(X, val, filepath_x, filepath_y, bag);
  X = readRowsFromBinary_umat(filepath_x, bag, 2);
  val = readRowsFromBinary(filepath_y, bag, 1);
  return List::create(X, val);
}

///
//write
///

// [[Rcpp::export]]
void writeBinaryVec(arma::vec data, std::string filepath) {
  std::ofstream out(filepath, std::ios::binary);
  if (!out.is_open()) {
    Rcpp::stop("Failed to open.");
  }
  out.write(reinterpret_cast<char*>(data.memptr()), data.n_elem * sizeof(double));
  out.close();
}

// [[Rcpp::export]]
void writeBinaryFile_umat(arma::umat data, std::string filepath) {
  std::ofstream out(filepath, std::ios::binary);
  if (!out.is_open()) {
    Rcpp::stop("Failed to open.");
  }
  out.write(reinterpret_cast<char*>(data.memptr()), data.n_elem * sizeof(unsigned int));
  out.close();
}