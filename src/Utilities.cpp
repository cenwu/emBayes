#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"Utilities.h"
#include<Rmath.h>
#include<iostream>
#include<stdio.h>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace R;

// [[Rcpp::export()]]
Rcpp::List logR(arma::vec& y, arma::mat& X, arma::mat& C, arma::vec alpha, arma::vec beta, double sigma2, double theta, double s0, double s1){

  arma::vec residule = y-X*beta-C*alpha;
  arma::vec dexps0 = exp(-abs(beta/s0))/(2*s0);
  arma::vec dexps1 = exp(-abs(beta/s1))/(2*s1);
  arma::vec numerator = theta*dexps1;
  arma::vec denominator = numerator + (1-theta)*dexps0 + 0.0000000001;
  arma::vec gamma = numerator/denominator;
  arma::vec invS = (1-gamma)/s0 + gamma/s1 + 0.0000000001;
  
  int n = X.n_rows;
  double p1 = 0.5*n*log(sigma2);
  double p2 = accu(residule%residule)/(2*sigma2);
  double p3 = accu(abs(beta)%invS);
  double p4 = accu(gamma*log(theta)+(1-gamma)*log(1-theta));
  double logver = -p1-p2-p3+p4-accu(pow(alpha,2)/2000)-log(sigma2);

  return Rcpp::List::create(Rcpp::Named("logver") = logver,
                            Rcpp::Named("Pgamma") = gamma,
                            Rcpp::Named("invS") = invS);
}
// [[Rcpp::export()]]
Rcpp::List logQR(arma::vec& y, arma::mat& X, arma::mat& C, arma::vec alpha, arma::vec beta, double sigma, double theta, double s0, double s1, double ep1, double ep22){
  int a = 1;
  int b = 1;
  int n = X.n_rows;
  arma::vec residule = y-X*beta-C*alpha;
  arma::vec dexps0 = exp(-abs(beta/s0))/(2*s0);
  arma::vec dexps1 = exp(-abs(beta/s1))/(2*s1);
  arma::vec numerator = theta*dexps1;
  arma::vec denominator = numerator + (1-theta)*dexps0 + 0.0000000001;
  arma::vec gamma = numerator/denominator;
  arma::vec invS = (1-gamma)/s0 + gamma/s1 + 0.0000000001;
  
  arma::vec delta2 = pow(residule,2)/(ep22*sigma);
  double gamma2 = (2+pow(ep1,2)/ep22)/sigma;
  
  arma::vec vn = y;
  arma::vec vp = y;
  
  for(int k = 0; k < n; k++){
    vn(k) = bessel_k(sqrt(gamma2*delta2(k)), 0.5-1 ,1)/bessel_k(sqrt(gamma2*delta2(k)), 0.5 ,1)*sqrt(gamma2/delta2(k));
    vp(k) = bessel_k(sqrt(gamma2*delta2(k)), 0.5+1 ,1)/bessel_k(sqrt(gamma2*delta2(k)), 0.5 ,1)*sqrt(delta2(k)/gamma2);
  }
    
  double p1 = 0.5*3*n*log(sigma);
  double p2 = accu((1/(2*sigma*ep22))*(vn%pow(residule,2)-2*residule*ep1+vp*pow(ep22,2)/4));
  double p3 = accu(abs(beta)%invS);
  double p4 = accu(gamma*log(theta)+(1-gamma)*log(1-theta));
  double logver = -p1-p2-p3+p4-accu(pow(alpha,2)/2000)+(a+1)*log(sigma)-(b/sigma);
  
  return Rcpp::List::create(Rcpp::Named("logver") = logver,
                            Rcpp::Named("Pgamma") = gamma,
                            Rcpp::Named("invS") = invS,
                            Rcpp::Named("vn") = vn,
                            Rcpp::Named("vp") = vp);
}



// [[Rcpp::export()]]
double Soft(double z, double lambda){
  if(z > lambda) return(z - lambda);
  else if(z < -lambda) return(z + lambda);
  else return(0);
}

