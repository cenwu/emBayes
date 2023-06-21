#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"EM.h"
#include"Utilities.h"
#include<Rmath.h>
#include<iostream>
#include<stdio.h>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace R;

// [[Rcpp::export()]]
Rcpp::List EMR(arma::vec& y, arma::mat& X, arma::mat& C, int n, int p, int q, arma::vec alpha, arma::vec beta, double sigma2, double theta, arma::vec Pgamma, arma::vec invS){
  
  theta = mean(Pgamma);
  
  arma::vec residule = y-X*beta-C*alpha;
  sigma2 = accu(residule%residule)/(n+2);
  
  for(int j = 0; j < p; j++){
    double denominator = accu(pow(X.col(j),2))/sigma2;
    double Zk = accu(X.col(j)%(residule+X.col(j)*beta(j)))/sigma2;
    beta(j) = Soft(Zk,invS(j))/denominator;
  }
  
  //alpha = accu(y-X*beta)/n;
  for(int k = 0; k < q; k++){
    double An = accu(C.col(k)%(residule+C.col(k)*alpha(k)));
    double Ad = accu(pow(C.col(k),2))+(2*sigma2/1000); 
    alpha(k) = An/Ad;
  }
  
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("theta") = theta);
}

// [[Rcpp::export()]]
Rcpp::List EMQR(arma::vec& y, arma::mat& X, arma::mat& C, int n, int p, int q, double quant, arma::vec alpha, arma::vec beta, double sigma, double theta, double s0, double s1, arma::vec Pgamma, arma::vec invS, double ep1, double ep22, arma::vec vn, arma::vec vp){
  
  theta = mean(Pgamma);
  int a = 1;
  int b = 1;
  arma::vec residule = y-X*beta-C*alpha;
  sigma = (accu(vn%pow(residule,2)-2*residule*ep1+vp*(pow(ep1,2)+2*ep22))+b)/((3*n+2*a+2)*ep22);
  
  for(int j = 0; j < p; j++){
    double d1 = ep22*sigma;
    double denominator = accu(vn%pow(X.col(j),2))/d1;
    
    double Z1 = accu(X.col(j)%(vn%(residule+X.col(j)*beta(j))))/d1;
    //double Z2 = accu(ep1*X.col(j))/ep22;
    double Z2 = accu(ep1*X.col(j))/d1;
    double Zk = Z1-Z2;
    beta(j) = Soft(Zk,invS(j))/denominator;
  }
  
  for(int k = 0; k < q; k++){
    double d2 = accu(vn%pow(C.col(k),2))+(ep22*sigma/1000);
    double A1 = accu(C.col(k)%(vn%(residule+C.col(k)*alpha(k))));
    double A2 = accu(ep1*C.col(k));    
    double Ak = A1 - A2;
    alpha(k) = Ak/d2;
  }
  //double a1 = accu(vn%(y-X*beta)-sigma*ep1);
  //double a2 = accu(vn);
  //alpha = a1/a2;
  
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("sigma") = sigma,
                            Rcpp::Named("theta") = theta);
}

