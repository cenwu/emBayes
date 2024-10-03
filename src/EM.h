#ifndef EM_h
#define EM_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List EMR(arma::vec& y, arma::mat& X, arma::mat& C, int n, int p, int q, arma::vec alpha, arma::vec beta, double sigma2, double theta, arma::vec Pgamma, arma::vec invS);
Rcpp::List EMQR(arma::vec& y, arma::mat& X, arma::mat& C, int n, int p, int q, double quant, arma::vec alpha, arma::vec beta, double sigma, double theta, double s0, double s1, arma::vec Pgamma, arma::vec invS, double ep1, double ep22, arma::vec vn, arma::vec vp);
Rcpp::List EMVCR(arma::vec& y, arma::mat& X, arma::mat& C, arma::mat& W, arma::mat& r, arma::mat& rs, int n, int m, arma::vec nt, int ns, int p, int q, arma::vec alpha, arma::vec beta, arma::vec g, double gn, double phi2, double sigma2, double theta, arma::vec Pgamma, arma::vec S);
Rcpp::List EMVCQR(arma::vec& y, arma::mat& X, arma::mat& C, arma::mat& W, arma::mat& r, arma::mat& rs, int n, int m, arma::vec nt, int ns, int p, int q, double quant, arma::vec alpha, arma::vec beta, arma::vec g, double gn, double phi2, double sigma, double theta, arma::vec Pgamma, arma::vec S, double ep1, double ep22, arma::vec vn, arma::vec vp);
#endif
