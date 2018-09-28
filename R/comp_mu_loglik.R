#' Calculate the Log-Likelihood of the COM-Poisson model
#' 
#' A function to compute the log-likelihood of the COM-Poisson model. 
#' 
#' @param param numeric vector: the model coefficients & the current value of \code{nu}. 
#' It is assumed that \code{nu} is in the last position of \code{param}.
#' @param y numeric vector: response variable
#' @param xx numeric matrix: the explanatory variables
#' @param offset numeric vector: a vector of lenght equal to the number of cases
#' 
#' @return 
#' The log-likelihood value of the COM-Poisson model.
comp_mu_loglik <-function(param, y, xx, offset){
  # compute negative loglikelihood for COMP-mu regression models
  # y is a n*1 column vector
  # xx is a n*q design matrix, including intercept
  # offset is a column vector matching the lenght of y
  n <- length(y)
  q <- ncol(xx)
  # regression coefficients
  #beta = param(1:q)  # not needed in loglikelihood
  lambda <- param[(q+1):(q+n)]
  # dispersion parameter
  nu <- param[(q+n+1)]
  # precompute quantities used later
  logfactorialy <- lgamma(y+1)
  Zcall <- Z(lambda, nu)
  meanlogfactorialy <- comp_mean_logfactorialy(lambda, nu)
  # compute loglikelihood 
  loglik <- sum(y*log(lambda) - nu*logfactorialy - log(Zcall))
  # compute gradient of negative loglikelihood (not currently used)
  #gradl <- -c(rep(0,q), y/lambda-comp_means(lambda,nu)/lambda,
  #            sum(-logfactorialy+meanlogfactorialy))
  return(loglik)
}
