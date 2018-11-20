// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;
using Eigen::MatrixXd;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
SEXP eigenMatMultipy(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenTCrossProd(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * t(B);
  
  return Rcpp::wrap(C);
}

// // [[Rcpp::export]]
// SEXP eigenSVD(const Eigen::Map<Eigen::MatrixXd> A_norm, int k, int q){
//   Eigen::MatrixXd C = A * t(B);
//   
//   return Rcpp::wrap(C);
// }
