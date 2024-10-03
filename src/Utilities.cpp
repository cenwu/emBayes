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
  double p2 = accu((1/(2*sigma*ep22))*(vn%pow(residule,2)-2*residule*ep1+vp*pow(ep1,2)))+accu(vp)/sigma-accu(log(vp))/2;
  double p3 = accu(abs(beta)%invS);
  double p4 = accu(gamma*log(theta)+(1-gamma)*log(1-theta));
  double logver = -p1-p2-p3+p4-accu(pow(alpha,2)/2000)-(a+1)*log(sigma)-(b/sigma);
  
  return Rcpp::List::create(Rcpp::Named("logver") = logver,
                            Rcpp::Named("Pgamma") = gamma,
                            Rcpp::Named("invS") = invS,
                            Rcpp::Named("vn") = vn,
                            Rcpp::Named("vp") = vp);
}

// [[Rcpp::export()]]
Rcpp::List logVCR(arma::vec& y, arma::mat& X, arma::mat& C, arma::mat& W, arma::mat& r, arma::mat& rs, int n, int m, int ns, arma::vec alpha, arma::vec beta, arma::vec g, double gn, double phi2, double sigma2, double theta, double s0, double s1){
  
  //int n = X.n_rows;
  //int m = W.n_cols;
  double a0 = 0.01;
  double b0 = 0.01;
  
  arma::vec dexps0(gn);
  arma::vec dexps1(gn);
  arma::vec normbeta(gn);
  
  for(int i = 0; i < gn; i++){
    double sub1 = accu(g.subvec(0,i));
    double sub2 = accu(g.subvec(0,i+1))-1;
    arma::vec betag = beta.subvec(sub1,sub2);
    normbeta(i) = accu(pow(betag,2));
    dexps0(i) = pow(s0,g(i+1))*exp(-s0*normbeta(i));
    dexps1(i) = pow(s1,g(i+1))*exp(-s1*normbeta(i));
  }
  arma::vec residule = y-X*beta-C*alpha-diagvec(W*r);
  arma::vec numerator = theta*dexps1;
  arma::vec denominator = numerator + (1-theta)*dexps0 + 0.0000000001;
  arma::vec gamma = numerator/denominator;
  arma::vec S = (1-gamma)*s0 + gamma*s1;
  
  double p1 = 0.5*n*log(sigma2);
  double p2 = accu(residule%residule)/(2*sigma2);
  double p3 = accu(S%normbeta);
  double p4 = accu(gamma*log(theta)+(1-gamma)*log(1-theta));
  double p5 = -m*ns*log(phi2)/2-accu(pow(rs,2))/(2*phi2)-(a0+1)*log(phi2)-b0/phi2;
  double logver = -p1-p2-p3+p4-accu(pow(alpha,2)/2000)-log(sigma2)+p5;
  
  return Rcpp::List::create(Rcpp::Named("r") = p5,
                            Rcpp::Named("logver") = logver,
                            Rcpp::Named("Pgamma") = gamma,
                            Rcpp::Named("S") = S);
}

// [[Rcpp::export()]]
Rcpp::List logVCQR(arma::vec& y, arma::mat& X, arma::mat& C, arma::mat& W, arma::mat& r, arma::mat& rs, int n, int m, int ns, arma::vec alpha, arma::vec beta, arma::vec g, double gn, double phi2, double sigma, double theta, double s0, double s1, double ep1, double ep22){
  int a = 1;
  int b = 1;
  //int n = X.n_rows;
  //int m = W.n_cols;
  double a0 = 0.01;
  double b0 = 0.01;
  
  arma::vec dexps0(gn);
  arma::vec dexps1(gn);
  arma::vec normbeta(gn);
  
  for(int i = 0; i < gn; i++){
    double sub1 = accu(g.subvec(0,i));
    double sub2 = accu(g.subvec(0,i+1))-1;
    arma::vec betag = beta.subvec(sub1,sub2);
    normbeta(i) = accu(pow(betag,2));
    dexps0(i) = pow(s0,g(i+1))*exp(-s0*normbeta(i));
    dexps1(i) = pow(s1,g(i+1))*exp(-s1*normbeta(i));
  }
  arma::vec residule = y-X*beta-C*alpha-diagvec(W*r);
  arma::vec numerator = theta*dexps1;
  arma::vec denominator = numerator + (1-theta)*dexps0 + 0.0000000001;
  arma::vec gamma = numerator/denominator;
  arma::vec S = (1-gamma)*s0 + gamma*s1;
  
  arma::vec delta2 = pow(residule,2)/(ep22*sigma);
  double gamma2 = (2+pow(ep1,2)/ep22)/sigma;
  
  arma::vec vn = y;
  arma::vec vp = y;
  
  for(int k = 0; k < n; k++){
    vn(k) = bessel_k(sqrt(gamma2*delta2(k)), 0.5-1 ,1)/bessel_k(sqrt(gamma2*delta2(k)), 0.5 ,1)*sqrt(gamma2/delta2(k));
    vp(k) = bessel_k(sqrt(gamma2*delta2(k)), 0.5+1 ,1)/bessel_k(sqrt(gamma2*delta2(k)), 0.5 ,1)*sqrt(delta2(k)/gamma2);
  }
  
  double p1 = 0.5*3*n*log(sigma);
  double p2 = accu((1/(2*sigma*ep22))*(vn%pow(residule,2)-2*residule*ep1+vp*pow(ep1,2)))+accu(vp)/sigma-accu(log(vp))/2;
  double p3 = accu(S%normbeta);
  double p4 = accu(gamma*log(theta)+(1-gamma)*log(1-theta));
  double p5 = -m*ns*log(phi2)/2-accu(pow(rs,2))/(2*phi2)-(a0+1)*log(phi2)-b0/phi2;
  double logver = -p1-p2-p3+p4-accu(pow(alpha,2)/2000)-(a+1)*log(sigma)-(b/sigma)+p5;
  
  return Rcpp::List::create(Rcpp::Named("r") = p3,
                            Rcpp::Named("logver") = logver,
                            Rcpp::Named("Pgamma") = gamma,
                            Rcpp::Named("S") = S,
                            Rcpp::Named("vn") = vn,
                            Rcpp::Named("vp") = vp);
} 



// [[Rcpp::export()]]
double Soft(double z, double lambda){
  if(z > lambda) return(z - lambda);
  else if(z < -lambda) return(z + lambda);
  else return(0);
}

