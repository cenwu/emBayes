#' @importFrom Rcpp sourceCpp
NULL

#' fit a model with given tuning parameters
#' 
#' This function performs penalized variable selection based on spike-and-slab quantile LASSO (ssQLASSO), spike-and-slab LASSO (ssLASSO), spike-and-slab quantile group LASSO varying coefficient mixed model (ssQVCM) and spike-and-slab group LASSO varying coefficient mixed model (ssVCM).
#' Typical usage is to first obtain the optimal spike scale and slab scale using cross-validation, then specify them in the 'emBayes' function.
#' @importFrom stats coefficients lm runif coef
#' @importFrom glmnet cv.glmnet glmnet
#' @param y a vector of response variable.
#' @param clin a matrix of clinical factors. It has default value NULL.
#' @param X a matrix of genetic factors.
#' @param W a matrix of random factors.
#' @param nt a vector of number of repeated measurements for each subject. They can be same or different.
#' @param group a vector of group sizes. They can be same or different.
#' @param quant value of quantile.
#' @param s0 value of the spike scale \eqn{s_{0}}.
#' @param s1 value of the slab scale \eqn{s_{1}}.
#' @param func methods to perform variable selection. Four choices are available. For non longitudinal analysis: "ssLASSO" and "ssQLASSO". For longitudinal varying-coefficient analysis: "ssVCM" and "ssQVCM".
#' @param error cutoff value for determining convergence. The algorithm reaches convergence if the difference in the expected log-likelihood of two iterations is less than the value of error. The default value is 0.01.
#' @param maxiter the maximum number of iterations that is used in the estimation algorithm. The default value is 200.
#' @details 
#' The current version of emBayes supports four types of methods: "ssLASSO", "ssQLASSO", "ssVCM" and "ssQVCM".
#' \itemize{
#' \item \strong{ssLASSO:} spike-and-slab LASSO fits a Bayesian linear regression through the EM algorithm.  
#' \item \strong{ssQLASSO:} spike-and-slab quantile LASSO fits a Bayesian quantile regression (based on asymmetric Laplace distribution) through the EM algorithm.
#' \item \strong{ssVCM:} spike-and-slab group LASSO varying coefficient mixed model fits a Bayesian linear mixed model through the EM algorithm.  
#' \item \strong{ssQVCM:} spike-and-slab quantile group LASSO varying coefficient mixed model fits a Bayesian quantile mixed model through the EM algorithm.
#' }
#' Users can choose the desired method by setting func="ssLASSO", "ssQLASSO", "ssVCM" or "ssQVCM".
#' @return A list with components:
#' \item{alpha}{a vector containing the estimated intercept and clinical coefficients.}
#' \item{intercept}{value of the estimated intercept.}
#' \item{clin.coe}{a vector of estimated clinical coefficients.}
#' \item{r}{a vector of estimated ranodm effect coefficients.}
#' \item{beta}{a vector of estimated beta coefficients.}
#' \item{sigma}{value of estimated asymmetric Laplace distribution scale parameter \eqn{\sigma}.}
#' \item{theta}{value of estimated probability parameter \eqn{\theta}.}
#' \item{iter}{value of number of iterations.}
#' \item{ll}{a vector of expectation of likelihood at each iteration.}
#' 
#' @examples 
#' data(data)
#' ##load the clinical factors, genetic factors, response and quantile data
#' clin=data$clin
#' X=data$X
#' y=data$y
#' quant=data$quant
#' 
#' ##generate tuning vectors of desired range 
#' t0 <- seq(0.01,0.015,length.out=2)
#' t1 <- seq(0.1,0.5,length.out=2)
#' 
#' ##perform cross-validation and obtain tuning parameters based on check loss
#' CV <- cv.emBayes(y,clin,X,W=NULL,nt=NULL,group=NULL,quant,t0,t1,k=5,
#' func="ssQLASSO",error=0.01,maxiter=200)
#' s0 <- CV$CL.s0
#' s1 <- CV$CL.s1
#' 
#' ##perform BQLSS under optimal tuning and calculate value of TP and FP for selecting beta 
#' EM <- emBayes(y,clin,X,W=NULL,nt=NULL,group=NULL,quant,s0,s1,func="ssQLASSO",
#' error=0.01,maxiter=200)
#' fit <- EM$beta
#' coef <- data$coef
#' tp <- sum(fit[coef!=0]!=0)
#' fp <- sum(fit[coef==0]!=0)
#' list(tp=tp,fp=fp)
#' 
#' @export

emBayes <- function(y,clin=NULL,X,W=NULL,nt=NULL,group=NULL,quant,s0,s1,func,error=0.01,maxiter=100){
  
  p <- ncol(X)
  n <- nrow(X)
  inter <- rep(1,n)
  #C <- cbind(inter,clin)
  #q <- ncol(C)
  ############################
  # func = "VC..."
  ns <- length(nt)
  nt <- c(0,nt)
  gn <- length(group)
  group <- c(0,group)
  if(is.null(W) == 0){
    W <- cbind(inter,W)
    m <- ncol(W)
  }  else{
    m=1
    W <- as.matrix(inter)
  }
  r <- matrix(0.1,m,n)
  rs <- matrix(0.1,m,ns)
  phi2=0.1
  ###########################
  
  #initial
  sigma=5
  sigma2=5
  theta=runif(1,0,1)
  lambda <- (cv.glmnet(X,y)$lambda.min)*1
  fit <- glmnet(X,y,lambda=lambda)
  beta <- coef(fit)[-1]
  
  y0 <- y-X%*%beta
  if(is.null(clin) == 0){
    C <- cbind(inter,clin)
    q <- ncol(C)
    regalpha <- lm(y0 ~ clin)
    alpha <- as.numeric(coefficients(regalpha))
  }
  else{
    C <- as.matrix(inter)
    q <- 1
    alpha <- mean(y0)
  }
  
  
  if(func=="ssQLASSO"){
    ep22 <- (2/(quant*(1-quant)))
    ep1 <- (1-2*quant)/(quant*(1-quant))
    
    ll <- c()
    loglike <- logQR(y,X,C,alpha,beta,sigma,theta,s0,s1,ep1,ep22)
    lk <- loglike$logver
    vn <- loglike$vn
    vp <- loglike$vp
    Pgamma <- loglike$Pgamma
    invS <- loglike$invS
    ll <- c(ll,lk)
    
    iter <- 0
    diff <- 1
    
    while( diff > error & iter < maxiter){
      iter  <- iter+1
      EM <- EMQR(y,X,C,n,p,q,quant,alpha,beta,sigma,theta,s0,s1,Pgamma,invS,ep1,ep22,vn,vp)
      alpha <- EM$alpha
      beta <- EM$beta
      sigma <- EM$sigma
      theta <- EM$theta
      
      loglike2 <- logQR(y,X,C,alpha,beta,sigma,theta,s0,s1,ep1,ep22)
      lk2 <- loglike2$logver
      ll <- c(ll,lk2)
      diff <- abs(lk2-lk)
      vn <- loglike2$vn
      vp <- loglike2$vp
      Pgamma <- loglike2$Pgamma
      invS <- loglike2$invS
      
      lk <- lk2
    }
  }
  
  else if(func=="ssLASSO"){
    ll <- c()
    loglike <- logR(y,X,C,alpha,beta,sigma2,theta,s0,s1)
    lk <- loglike$logver
    ll <- c(ll,lk)
    Pgamma <- loglike$Pgamma
    invS <- loglike$invS
    
    iter <- 0
    diff <- 1
    
    while( diff > error & iter < maxiter){
      iter  <- iter+1
      EM <- EMR(y,X,C,n,p,q,alpha,beta,sigma2,theta,Pgamma,invS)
      alpha <- EM$alpha
      beta <- EM$beta
      sigma2 <- EM$sigma2
      theta <- EM$theta
      
      loglike2 <- logR(y,X,C,alpha,beta,sigma2,theta,s0,s1)
      lk2 <- loglike2$logver
      ll <- c(ll,lk2)
      diff <- abs(lk2-lk)
      Pgamma <- loglike2$Pgamma
      invS <- loglike2$invS
      
      lk <- lk2
    }
  }
  
  else if(func=="ssVCM"){
    ll <- c()
    loglike <- logVCR(y,X,C,W,r,rs,n,m,ns,alpha,beta,group,gn,phi2,sigma2,theta,s0,s1)
    lk <- loglike$logver
    ll <- c(ll,lk)
    Pgamma <- loglike$Pgamma
    S <- loglike$S
    
    iter <- 0
    diff <- 1
    d <- c()
    
    while( diff > error & iter < maxiter){
      iter  <- iter+1
      EM <- EMVCR(y,X,C,W,r,rs,n,m,nt,ns,p,q,alpha,beta,group,gn,phi2,sigma2,theta,Pgamma,S)
      alpha <- EM$alpha
      beta <- EM$beta
      r <- EM$r
      
      for(i in 1:ns){
        t1 <- sum(nt[1:i])+1
        t2 <- sum(nt[1:(i+1)])
        rs[,i] <- unique(r[,t1:t2],MARGIN=2)
      }
      
      phi2 <- EM$phi2
      sigma2 <- EM$sigma2
      theta <- EM$theta
      
      loglike2 <- logVCR(y,X,C,W,r,rs,n,m,ns,alpha,beta,group,gn,phi2,sigma2,theta,s0,s1)
      lk2 <- loglike2$logver
      ll <- c(ll,lk2)
      diff <- abs(lk2-lk)
      d <- c(d,diff)
      #if(is.null(diff) == 0){break}
      Pgamma <- loglike2$Pgamma
      S <- loglike2$S
      
      lk <- lk2
    }
  }
  
  else{
    ll <- c()
    loglike <- logVCQR(y,X,C,W,r,rs,n,m,ns,alpha,beta,group,gn,phi2,sigma,theta,s0,s1,ep1,ep22)
    lk <- loglike$logver
    vn <- loglike$vn
    vp <- loglike$vp
    Pgamma <- loglike$Pgamma
    S <- loglike$S
    ll <- c(ll,lk)
    
    iter <- 0
    diff <- 1
    d <- c(diff)
    ph <- c()
    
    while( diff > error  & iter < maxiter){
      iter  <- iter+1
      cat(iter)
      EM <- EMVCQR(y,X,C,W,r,rs,n,m,nt,ns,p,q,quant,alpha,beta,group,gn,phi2,sigma,theta,Pgamma,S,ep1,ep22,vn,vp)
      alpha <- EM$alpha
      beta2 <- EM$beta
      r <- EM$r
      
      for(i in 1:ns){
        t1 <- sum(nt[1:i])+1
        t2 <- sum(nt[1:(i+1)])
        rs[,i] <- unique(r[,t1:t2],MARGIN=2)
      }
      
      phi2 <- EM$phi2
      ph <- c(ph,phi2)
      sigma <- EM$sigma
      theta <- EM$theta
      
      loglike2 <- logVCQR(y,X,C,W,r,rs,n,m,ns,alpha,beta,group,gn,phi2,sigma,theta,s0,s1,ep1,ep22)
      lk2 <- loglike2$logver
      ll <- c(ll,lk2)
      #diff <- abs(lk2-lk)
      diff <- sum(abs(beta-beta2))
      d <- c(d,diff)
      #if(is.null(diff) == 0){break}
      vn <- loglike2$vn
      vp <- loglike2$vp
      Pgamma <- loglike2$Pgamma
      S <- loglike2$S
      
      lk <- lk2
      beta <- beta2
    }
  }
    
  intercept <- alpha[1]
  clin.coe <- alpha[-1]
  this.call = match.call()

  return(list(call = this.call,alpha=alpha,r=r,intercept=intercept,clin.coe=clin.coe,beta=beta,sigma=sigma,theta=theta,iter=iter,ll=ll))
}