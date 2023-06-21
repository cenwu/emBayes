#ifndef Utilities_h
#define Utilities_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List logR(arma::vec& y, arma::mat& X, arma::mat& C, arma::vec alpha, arma::vec beta, double sigma2, double theta, double s0, double s1);
Rcpp::List logQR(arma::vec& y, arma::mat& X, arma::mat& C, arma::vec alpha, arma::vec beta, double sigma, double theta, double s0, double s1, double ep1, double ep22);
double Soft(double z, double lambda);

#endif
