
<!-- README.md is generated from README.Rmd. Please edit that file -->

# emBayes

Variable selection methods have been extensively developed for analyzing high-dimensional omics data within both the frequentist and Bayesian frameworks. This package implemented the spike-and-slab quantile LASSO which has been developed along the line of Bayesian hierarchical model but deeply rooted in the frequentist regularization methods by utilizing the Expectation–Maximization (EM) algorithm. Therefore, the proposed method borrows strength from both the frequentist and Bayesian frameworks while overcoming their respective limitations. The spike-and-slab quantile LASSO can handle data irregularity in terms of skewness and outliers in the disease trait, compared to its nonrobust alternative, the spike-and-slab LASSO, which has also been implemented in the package. The core module of this package is developed in C++. 

## How to install

    - To install from Github, run these two lines of code in R

<!-- end list -->

    install.packages("devtools")
    devtools::install_github("cenwu/emBayes")

  - Released versions of mixedBayes are available on CRAN
    [(link)](https://cran.r-project.org/package=emBayes), and can be
    installed within R via

<!-- end list -->

    install.packages("mixedBayes")

## Example

    library(emBayes)
    data(genes)
    ## Load the clinical factors, genetic factors, response and quantile data
    clin=genes$clin
    X=genes$X
    y=genes$y
    quant=genes$quant
    ## Generate tuning vectors of desired range
    t0 <- seq(0.01,0.015,length.out=2)
    t1 <- seq(0.1,0.5,length.out=2)
    ## Perform cross-validation and obtain tuning parameters based on check loss
    CV <- cv.emBayes(y,clin,X,quant,t0,t1,k=5,func="BQLSS",error=0.01,maxiter=2000)
    s0 <- CV$CL.s0
    s1 <- CV$CL.s1
    ## Perform BQLSS under optimal tuning and calculate value of TP and FP for selecting beta
    EM <- emBayes(y,clin,X,quant,s0,s1,func="BQLSS",error=0.01,maxiter=2000)
    fit <- EM$beta
    coef <- genes$coef
    tp <- sum(fit[coef!=0]!=0)
    fp <- sum(fit[coef==0]!=0)
    list(tp=tp,fp=fp)



## Methods

This package provides implementation for methods proposed in

   - Liu, Y., Ren, J., Ma, S., Wang, W. and Wu, C. (2023+). The Spike-and-Slab Quantile LASSO.(submitted)

