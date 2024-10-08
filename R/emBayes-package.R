#' @useDynLib emBayes, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @docType package
#' @keywords overview
#' @name emBayes-package
#' @title Robust Bayesian Variable Selection via Expectation-Maximization
#' @aliases emBayes-package
#' @description
#' This package provides the implementation of the spike-and-slab quantile LASSO (ssQLASSO) and spike-and-slab quantile group LASSO varying coefficient mixed model (ssQVCM) which combines the strength of Bayesian robust variable selection and the Expectation-Maximization (EM) coordinate descent approach. The alternative methods spike-and-slab LASSO (ssLASSO) and spike-and-slab group LASSO varying coefficient mixed model (ssVCM) are also included in the package.
#' 
#' @details
#' Two user friendly, integrated interface \strong{cv.emBayes()} and \strong{emBayes()} allows users to flexibly choose the variable selection method by specifying the following parameter:
#' \tabular{rl}{
#' quant: \tab to specify different quantiles when using robust methods.\cr\cr
#' func: \tab the model to perform variable selection. Four choices are available: \cr \tab "ssLASSO", "ssQLASSO", "ssVCM" and "ssQVCM". \cr\cr
#' error: \tab to specify the difference between expectations of likelihood of two  \cr \tab consecutive iterations. It can be used to determine convergence.\cr\cr
#' maxiter: \tab to specify the maximum number of iterations.
#' }
#' 
#' Function cv.emBayes() returns cross-validation errors based on the check loss, least squares loss and Schwarz Information Criterion along with the corresponding optimal tuning parameters.
#' Function emBayes() returns the estimated intercept, clinical coefficients, beta coefficients, scale parameter, probability parameter, number of iterations and expectation of likelihood at each iteration. 
#' 
#' @references
#' Liu, Y., Ren, J., Ma, S., and Wu, C. (2024). The Spike-and-Slab Quantile LASSO for Robust Variable Selection in Cancer Genomics Studies.
#' \emph{Statistics in Medicine.}
#' 
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y., and Wu, C. (2022). Robust Bayesian variable selection for gene–environment interactions.
#' \emph{Biometrics.} \doi{10.1111/biom.13670}
#' 
#' Ren, J., Du, Y., Li, S., Ma, S., Jiang,Y. and Wu, C. (2019). Robust network-based regularization
#' and variable selection for high dimensional genomics data in cancer prognosis.
#' {\emph{Genet. Epidemiol.}, 43:276-291} \doi{10.1002/gepi.22194}
#' 
#' Wu, C., Zhang, Q., Jiang,Y. and Ma, S. (2018). Robust network-based analysis of the associations between (epi)genetic measurements.
#' {\emph{J Multivar Anal.}, 168:119-130} \doi{10.1016/j.jmva.2018.06.009}
#' 
#' Tang, Z., Shen, Y., Zhang, X., and Yi, N. (2017). The spike-and-slab lasso generalized linear models for prediction and associated genes detection. 
#' {\emph{Genetics}, 205(1), 77-88} \doi{10.1534/genetics.116.192195}
#' 
#' Tang, Z., Shen, Y., Zhang, X., and Yi, N. (2017). The spike-and-slab lasso Cox model for survival prediction and associated genes detection. 
#' {\emph{Bioinformatics}, 33(18), 2799-2807} \doi{10.1093/bioinformatics/btx300}
#' 
#' Wu, C., and Ma, S. (2015). A selective review of robust variable selection with applications in bioinformatics.
#' {\emph{Briefings in Bioinformatics}, 16(5), 873–883} \doi{10.1093/bib/bbu046} 
#' 
#' Zhou, Y. H., Ni, Z. X., and Li, Y. (2014). Quantile regression via the EM algorithm. 
#' {\emph{Communications in Statistics-Simulation and Computation}, 43(10), 2162-2172} \doi{10.1080/03610918.2012.746980}
#' 
#' Ročková, V., and George, E. I. (2014). EMVS: The EM approach to Bayesian variable selection. 
#' {\emph{Journal of the American Statistical Association}, 109(506), 828-846} \doi{10.1080/01621459.2013.869223}
#' 
#' Li, Q., Lin, N., and Xi, R. (2010). Bayesian regularized quantile regression.
#' {\emph{Bayesian Analysis}, 5(3), 533-556} \doi{10.1214/10-BA521}
#' 
#' George, E. I., and McCulloch, R. E. (1993). Variable selection via Gibbs sampling. 
#' {\emph{Journal of the American Statistical Association}, 88(423), 881-889} \doi{10.1080/01621459.1993.10476353}
#' 
#' @seealso \code{\link{cv.emBayes}} \code{\link{emBayes}}
NULL

