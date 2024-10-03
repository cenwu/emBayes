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
    arma::vec residule1 = y-X*beta-C*alpha;
    double denominator = accu(pow(X.col(j),2))/sigma2;
    double Zk = accu(X.col(j)%(residule1+X.col(j)*beta(j)))/sigma2;
    beta(j) = Soft(Zk,invS(j))/denominator;
  }
  
  //alpha = accu(y-X*beta)/n;
  for(int k = 0; k < q; k++){
    arma::vec residule2 = y-X*beta-C*alpha;
    double An = accu(C.col(k)%(residule2+C.col(k)*alpha(k)));
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
    arma::vec residule1 = y-X*beta-C*alpha;
    double d1 = ep22*sigma;
    double denominator = accu(vn%pow(X.col(j),2))/d1;
    
    double Z1 = accu(X.col(j)%(vn%(residule1+X.col(j)*beta(j))))/d1;
    //double Z2 = accu(ep1*X.col(j))/ep22;
    double Z2 = accu(ep1*X.col(j))/d1;
    double Zk = Z1-Z2;
    beta(j) = Soft(Zk,invS(j))/denominator;
  }
  
  for(int k = 0; k < q; k++){
    arma::vec residule2 = y-X*beta-C*alpha;
    double d2 = accu(vn%pow(C.col(k),2))+(ep22*sigma/1000);
    double A1 = accu(C.col(k)%(vn%(residule2+C.col(k)*alpha(k))));
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


// [[Rcpp::export()]]
Rcpp::List EMVCR(arma::vec& y, arma::mat& X, arma::mat& C, arma::mat& W, arma::mat& r, arma::mat& rs, int n, int m, arma::vec nt, int ns, int p, int q, arma::vec alpha, arma::vec beta, arma::vec g, double gn, double phi2, double sigma2, double theta, arma::vec Pgamma, arma::vec S){
  
  double a0 = 0.01;
  double b0 = 0.01;
  
  theta = mean(Pgamma);
  
  arma::vec residule = y-X*beta-C*alpha-diagvec(W*r);
  sigma2 = accu(residule%residule)/(n+2);
  
  double ph1 = accu(pow(rs,2))/2+b0;
  double ph2 = ns*m/2+a0+1;
  phi2 = ph1/ph2;
  
  for(int i = 0; i < gn; i++){
    double sub1 = accu(g.subvec(0,i));
    double sub2 = accu(g.subvec(0,i+1))-1;
    arma::vec betag = beta.subvec(sub1,sub2);
    arma::mat Xg = X.cols(sub1,sub2);
    arma::vec rg = y-X*beta-C*alpha + Xg*betag-diagvec(W*r);
    arma::vec ug = inv(Xg.t()*Xg)*Xg.t()*rg;
    
    arma::vec tg = Xg*ug;
    double norm = accu(pow(tg,2));
    
    beta.subvec(sub1,sub2) = Soft(norm,sigma2*S(i))/norm*ug;
  }
  
  //int ns = n/nt;
  for(int j = 0; j < ns; j++){
    double t1 = accu(nt.subvec(0,j));
    double t2 = accu(nt.subvec(0,j+1))-1;
    //int t1 = j*nt;
    //int t2 = (j+1)*nt-1;
    arma::mat Wj = W.rows(t1,t2);
    arma::mat rj = r.cols(t1,t2);
    arma::vec res = y.subvec(t1,t2)-X.rows(t1,t2)*beta-C.rows(t1,t2)*alpha;
    for(int l = 0; l < m; l++){
      arma::vec Wjl = Wj.col(l);
      double r1 = accu((res-diagvec(Wj*rj))%Wjl)/sigma2;
      double r2 = accu(pow(Wjl,2))/sigma2+(1/phi2);
      
      for(int t0 = 0; t0 < nt(j+1); t0++){
        ((r.cols(t1,t2)).row(l))(t0) = r1/r2;
      }
    }
  }
  
  for(int k = 0; k < q; k++){
    arma::vec residule2 = y-X*beta-C*alpha-diagvec(W*r);
    double An = accu(C.col(k)%(residule2+C.col(k)*alpha(k)));
    double Ad = accu(pow(C.col(k),2))+(2*sigma2/1000); 
    alpha(k) = An/Ad;
  }
  
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("r") = r,
                            Rcpp::Named("phi2") = phi2,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("theta") = theta);
}

// [[Rcpp::export()]]
Rcpp::List EMVCQR(arma::vec& y, arma::mat& X, arma::mat& C, arma::mat& W, arma::mat& r, arma::mat& rs, int n, int m, arma::vec nt, int ns, int p, int q, double quant, arma::vec alpha, arma::vec beta, arma::vec g, double gn, double phi2, double sigma, double theta, arma::vec Pgamma, arma::vec S, double ep1, double ep22, arma::vec vn, arma::vec vp){
  
  double a0 = 0.01;
  double b0 = 0.01;
  int a = 1;
  int b = 1;
  
  theta = mean(Pgamma);
  
  arma::vec residule = y-X*beta-C*alpha-diagvec(W*r);
  sigma = (accu(vn%pow(residule,2)-2*residule*ep1+vp*(pow(ep1,2)+2*ep22))+b)/((3*n+2*a+2)*ep22);
  
  //double ph1 = accu(pow(r,2))/(2*nt)+b0;
  //double ph2 = n*m/(nt*2)+a0+1;
  //phi2 = ph1/ph2;
  
  arma::mat dv = diagmat(vn);
  arma::vec o(n,arma::fill::ones);
  arma::vec oo(5,arma::fill::zeros);
  //o.subvec(0,4) = oo;
  for(int i = 0; i < gn; i++){
    double sub1 = accu(g.subvec(0,i));
    double sub2 = accu(g.subvec(0,i+1))-1;
    arma::vec betag = beta.subvec(sub1,sub2);
    arma::mat Xg = X.cols(sub1,sub2);
    arma::vec rg = y-X*beta-C*alpha+Xg*betag-diagvec(W*r);
    arma::vec ug = inv(Xg.t()*dv*Xg)*(Xg.t()*dv*rg-ep1*Xg.t()*o);
    
    arma::vec tg = sqrt(dv)*Xg*ug;
    double norm = accu(pow(tg,2));
    double t0 = ep22*sigma*S(i);
    
    beta.subvec(sub1,sub2) = Soft(norm,t0)/norm*ug;
  }
  
  //int ns = n/nt;
  for(int j = 0; j < ns; j++){
    double t1 = accu(nt.subvec(0,j));
    double t2 = accu(nt.subvec(0,j+1))-1;
    //int t1 = j*nt;
    //int t2 = (j+1)*nt-1;
    arma::vec res = y.subvec(t1,t2)-X.rows(t1,t2)*beta-C.rows(t1,t2)*alpha;
    arma::mat Wj = W.rows(t1,t2);
    arma::mat rj = r.cols(t1,t2);
    arma::vec vnj = vn.subvec(t1,t2);
    for(int l = 0; l < m; l++){
      arma::vec Wjl = Wj.col(l);
      double r1 = accu((res-diagvec(Wj*rj))%Wjl%vnj)/(sigma*ep22)-ep1*accu(Wjl)/(ep22*sigma);
      double r2 = accu(pow(Wjl,2)%vnj)/(sigma*ep22)+(1/phi2);
      for(int t0 = 0; t0 < nt(j+1); t0++){
        ((r.cols(t1,t2)).row(l))(t0) = r1/r2;
      }
    }
  }
  
  double ph1 = accu(pow(rs,2))/2+b0;
  double ph2 = ns*m/2+a0+1;
  phi2 = ph1/ph2;
  
  for(int k = 0; k < q; k++){
    arma::vec residule2 = y-X*beta-C*alpha-diagvec(W*r);
    double d2 = accu(vn%pow(C.col(k),2))+(ep22*sigma/1000);
    double A1 = accu(C.col(k)%(vn%(residule2+C.col(k)*alpha(k))));
    double A2 = accu(ep1*C.col(k));    
    double Ak = A1 - A2;
    alpha(k) = Ak/d2;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("r") = r,
    Rcpp::Named("phi2") = phi2,
    Rcpp::Named("sigma") = sigma,
    Rcpp::Named("theta") = theta);

}
