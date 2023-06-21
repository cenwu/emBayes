#ifndef EM_h
#define EM_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List EMR(arma::vec& y, arma::mat& X, arma::mat& C, int n, int p, int q, arma::vec alpha, arma::vec beta, double sigma2, double theta, arma::vec Pgamma, arma::vec invS);
Rcpp::List EMQR(arma::vec& y, arma::mat& X, arma::mat& C, int n, int p, int q, double quant, arma::vec alpha, arma::vec beta, double sigma, double theta, double s0, double s1, arma::vec Pgamma, arma::vec invS, double ep1, double ep22, arma::vec vn, arma::vec vp);
#endif
