#' simulated gene expression example data 
#' 
#' Simulated gene expression data for demonstrating the usage of emBayes.
#' 
#' @docType data
#' @keywords datasets
#' @name data
#' @usage data(data)
#' @format The data file consists of five components: y, clin, X, quant, coef and clin.coe. The coefficients and clinical coefficients are the true values of parameters used for generating response y. They can be used for performance evaluation.
#' 
#' @details
#'
#' \strong{The data model for generating response}
#'
#'Let \eqn{y_{i}} be the response of the \eqn{i}-th subject (1\eqn{\leq} i\eqn{\leq} n). We have \eqn{z_{i}=(1,z_{i1},\dots,z_{iq})^{\top}} being a \eqn{(q+1)}-dimensional vector of which the last \eqn{q} components indicate clinical factors and \eqn{x_{i}=(x_{i1},\dots,x_{ip})^{\top}} denoting a \eqn{p}-dimensional vector of genetic factors. The linear quantile regression model for the \eqn{\tau}-th quantile \eqn{(0<\tau<1)} is:
#'
#'\deqn{y_i=z_i^\top\alpha+x_i^\top\beta+\epsilon_i}
#'where \eqn{\alpha=(\alpha_0,\cdots,\alpha_q)^\top} contains the intercept and the regression coefficients for the clinical covariates. \eqn{\beta=(\beta_1,\cdots,\beta_p)^\top} are the regression coefficients and random error \eqn{\epsilon_{i}=(\epsilon_{1},...,\epsilon_{n})^\top} is set to follow a T2 distribution and has value \eqn{0} at its \eqn{\tau}-th quantile. 
#' @seealso \code{\link{emBayes}}
NULL
