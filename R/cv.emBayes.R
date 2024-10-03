#' @importFrom Rcpp sourceCpp
NULL

#' k-folds cross-validation for 'emBayes'
#' 
#' This function performs cross-validation and returns the optimal values of the tuning parameters.
#' @param y a vector of response variable.
#' @param clin a matrix of clinical factors. It has default value NULL.
#' @param X a matrix of genetic factors.
#' @param W a matrix of random factors.
#' @param nt a vector of number of repeated measurements for each subject. They can be same or different.
#' @param group a vector of group sizes. They can be same or different.
#' @param quant value of quantile.
#' @param t0 a user-supplied sequence of the spike scale \eqn{s_{0}}.
#' @param t1 a user-supplied sequence of the slab scale \eqn{s_{1}}.
#' @param k number of folds for cross-validation.
#' @param func methods to perform variable selection. Four choices are available. For non longitudinal analysis: "ssLASSO" and "ssQLASSO". For longitudinal varying-coefficient analysis: "ssVCM" and "ssQVCM".
#' @param error cutoff value for determining convergence. The algorithm reaches convergence if the difference in the expected log-likelihood of two iterations is less than the value of error. The default value is 0.01.
#' @param maxiter the maximum number of iterations that is used in the estimation algorithm. The default value is 200.
#' @details
#' When performing cross-validation for emBayes, function cv.emBayes returns two sets of optimal tuning parameters and their corresponding cross-validation error matrices. 
#' The spike scale parameter \eqn{CL.s0} and the slab scale parameter \eqn{CL.s1} are obtained based on the quantile check loss. 
#' The spike scale parameter \eqn{SL.s0} and the slab scale parameter \eqn{SL.s1} are obtained based on the least squares loss. 
#' The spike scale parameter \eqn{SIC.s0} and the slab scale parameter \eqn{SIC.s1} are obtained based on the Schwarz Information Criterion (SIC).
#' Corresponding error matrices \eqn{CL.CV}, \eqn{SL.CV} and \eqn{SIC.CV} can also be obtained from the output.
#' 
#' Schwarz Information Criterion has the following form:
#' \deqn{SIC=\log\sum_{i=1}^nL(y_i-\hat{y_i})+\frac{\log n}{2n}edf}
#' where \eqn{L(\cdot)} is the check loss and \eqn{edf} is the number of close to zero residuals \eqn{(\leq 0.001)}. For non-robust method “ssLASSO”, one should use least squares loss for tuning selection. For robust method “ssQLASSO”, one can either use quantile check loss or SIC for tuning selection. We suggest using SIC, since it has been extensively utilized for tuning selection in high-dimensional quantile regression, as documented in numerous literature sources.
#' 
#' @return A list with components:
#' \item{CL.s0}{the optimal spike scale under check loss.}
#' \item{CL.s1}{the optimal slab scale under check loss.}
#' \item{SL.s0}{the optimal slab scale under least squares loss.}
#' \item{SL.s1}{the optimal slab scale under least squares loss.}
#' \item{SIC.s0}{the optimal slab scale under SIC.}
#' \item{SIC.s1}{the optimal slab scale under SIC.}
#' \item{CL.CV}{cross-validation error matrix under check loss.}
#' \item{SL.CV}{cross-validation error matrix under least squares loss.}
#' \item{SIC.CV}{cross-validation error matrix under SIC.}
#' 
#' @export

cv.emBayes <- function(y,clin=NULL,X,W=NULL,nt=NULL,group=NULL,quant,t0,t1,k,func,error=0.01,maxiter=100){
  l1 <- length(t0)
  l2 <- length(t1)
  n <- nrow(X)
  inter <- rep(1,n)
  C <- cbind(inter,clin)
  
  CV1 <- matrix(0,l1,l2)
  CV2 <- matrix(0,l1,l2)
  CV3 <- matrix(0,l1,l2)
  
  if(grepl("VC",func)==0){
    s <- sample(1:n,n,replace=FALSE)
    folds <- cut(s,breaks=k,labels=FALSE)
  
    for(i in 1:l1){
        for(j in 1:l2){
          rob <- rep(0,k)
          nrob <- rep(0,k)
          SIC <- rep(0,k)
          for(f in 1:k){
            k.sub <- which(folds==f)
            EM <- emBayes(y[-k.sub], clin[-k.sub,], X[-k.sub,],W=NULL,nt=NULL,group=NULL,quant,t0[i],t1[j],func,error,maxiter)
            beta_k <- EM$beta
            alpha_k <- EM$alpha
            residual <- y[k.sub]-X[k.sub,]%*%beta_k-C[k.sub,]%*%alpha_k
            edf <- sum(abs(residual) <= 1e-03)
            rob[f] <- sum(residual*(quant-(residual < 0))) 
            nrob[f] <- sum(residual^2)
            SIC[f] <- log(rob[f])+log(n)*edf/(2*n)
          }
          CV1[i,j] <- sum(rob)/k
          CV2[i,j] <- sum(nrob)/k
          CV3[i,j] <- sum(SIC)/k
        }
    }
  }
  
  else{
    ns <- length(nt)
    s <- sample(1:ns,ns,replace=FALSE)
    folds <- cut(s,breaks=k,labels=FALSE)
    
    for(i in 1:l1){
      for(j in 1:l2){
        rob <- rep(0,k)
        nrob <- rep(0,k)
        SIC <- rep(0,k)
        for(f in 1:k){
          k.sub <- which(folds==f)
          EM <- emBayes(y[-k.sub], clin[-k.sub,], X[-k.sub,],W[-k.sub,],nt[-k.sub],group,quant,t0[i],t1[j],func,error,maxiter)
          beta_k <- EM$beta
          alpha_k <- EM$alpha
          r_k <- EM$r
          residual <- y[k.sub]-X[k.sub,]%*%beta_k-C[k.sub,]%*%alpha_k-diag(W[k.sub,]%*%r_k)
          edf <- sum(abs(residual) <= 1e-03)
          rob[f] <- sum(residual*(quant-(residual < 0))) 
          nrob[f] <- sum(residual^2)
          SIC[f] <- log(rob[f])+log(n)*edf/(2*n)
        }
        CV1[i,j] <- sum(rob)/k
        CV2[i,j] <- sum(nrob)/k
        CV3[i,j] <- sum(SIC)/k
      }
    }
  }
  

  indices1 <- which(CV1 == min(CV1), arr.ind=TRUE)
  indices2 <- which(CV2 == min(CV2), arr.ind=TRUE)
  indices3 <- which(CV3 == min(CV3), arr.ind=TRUE)
  return(list("CL.s0"=t0[indices1[1]],"CL.s1"=t1[indices1[2]],"CL.CV"=CV1,"SL.s0"=t0[indices2[1]],"SL.s1"=t1[indices2[2]],"SL.CV"=CV2,"SIC.s0"=t0[indices3[1]],"SIC.s1"=t1[indices3[2]],"SIC.CV"=CV3))
}