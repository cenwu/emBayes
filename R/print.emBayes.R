#' print an emBayes result
#'
#' Print a summary of an 'emBayes' result
#'
#' @param x emBayes result
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @return Print a list of output from a emBayes object.
#' @seealso \code{\link{emBayes}}
#' @export
print.emBayes=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nCall: ", deparse(x$call), "\n")
  cat("\nAlpha:\n")
  print(x$alpha, digits)
  cat("\nIntercept:\n")
  print(x$intercept, digits)
  cat("\nClinical Coefficients:\n")
  print(x$clin.coe, digits)
  cat("\nBeta Coefficients:\n")
  print(x$beta, digits)
  cat("\nSigma:\n")
  print(x$sigma, digits)
  cat("\nTheta:\n")
  print(x$theta, digits)
  cat("\nIterations:\n")
  print(x$iter, digits)
  cat("\nLog Likelihood:\n")
  print(x$ll, digits)
}



#' print an cv.emBayes result
#'
#' Print a summary of an 'cv.emBayes' result
#'
#' @param x cv.emBayes result
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @return Print a list of output from a cv.emBayes object.
#' @seealso \code{\link{cv.emBayes}}
#' @export
print.cv.emBayes=function(x, digits = max(3, getOption("digits") - 3),...){
  #cat("\nCall: ", deparse(x$call), "\n")
  cat("\ncheck loss spike tuning:\n")
  print(x$CL.s0, digits)
  cat("\ncheck loss slab tuning:\n")
  print(x$CL.s1, digits)
  cat("\nleast squares loss spike tuning:\n")
  print(x$SL.s0, digits)
  cat("\nleast squares loss slab tuning:\n")
  print(x$SL.s1, digits)
  cat("\nSIC spike tuning:\n")
  print(x$SIC.s0, digits)
  cat("\nSIC slab tuning:\n")
  print(x$SIC.s1, digits)
  cat("\ncheck loss CV error:\n")
  print(x$CL.CV, digits)
  cat("\nleast squares loss CV error:\n")
  print(x$SL.CV, digits)
  cat("\nSIC CV error:\n")
  print(x$SIC.CV, digits)
}