#ifndef Utilities_h
#define Utilities_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List logR(arma::vec& y, arma::mat& X, arma::mat& C, arma::vec alpha, arma::vec beta, double sigma2, double theta, double s0, double s1);
Rcpp::List logQR(arma::vec& y, arma::mat& X, arma::mat& C, arma::vec alpha, arma::vec beta, double sigma, double theta, double s0, double s1, double ep1, double ep22);
Rcpp::List logVCR(arma::vec& y, arma::mat& X, arma::mat& C, arma::mat& W, arma::mat& r, arma::mat& rs, int n, int m, int ns, arma::vec alpha, arma::vec beta, arma::vec g, double gn, double phi2, double sigma2, double theta, double s0, double s1);
Rcpp::List logVCQR(arma::vec& y, arma::mat& X, arma::mat& C, arma::mat& W, arma::mat& r, arma::mat& rs, int n, int m, int ns, arma::vec alpha, arma::vec beta, arma::vec g, double gn, double phi2, double sigma, double theta, double s0, double s1, double ep1, double ep22);
double Soft(double z, double lambda);

#endif
