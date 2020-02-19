#' Calculate the Log-Likelihood of the COM-Poisson model
#' 
#' A function to compute the log-likelihood of the COM-Poisson model. 
#' 
#' @param param numeric vector: the model coefficients & the current value of \code{nu}. 
#' It is assumed that \code{nu} is in the last position of \code{param}.
#' @param y numeric vector: response variable
#' @param xx numeric matrix: the explanatory variables
#' @param offset numeric vector: a vector of length equal to the number of cases
#' @param summax maximum number of terms to be considered in the truncated sum
#' @param log_nu numeric: nu in log-scale
#' @param mu numeric vector: fitted mean parameters
#' @return 
#' \code{comp_mu_loglik} returns the log-likelihood value of the COM-Poisson model based on Huang (2018).
#' \code{comp_mu_neg_loglik_log_nu_only} returns the negative log-likelihood value of the COM-Poisson model based on Ribeiro Jr et al. (2018)'s specification to use in conjunction with \code{optim}.
#' @name comp_mu_loglik
NULL

#' @rdname comp_mu_loglik
comp_mu_loglik <-function(param, y, xx, offset, summax){
  # compute loglikelihood for COMP-mu regression models
  # y is a n*1 column vector
  # xx is a n*q design matrix, including intercept
  # offset is a column vector matching the length of y
  n <- length(y)
  q <- ncol(xx)
  # regression coefficients
  beta = param[1:q] 
  eta <- t(xx%*%beta)[1,]
  mu <- exp(eta+offset)
  lambda <- param[(q+1):(q+n)]
  # dispersion parameter
  if (length(param)==q+n+1){
    nu <- param[(q+n+1)]} else {
      nu <- param[(q+n+1):(q+2*n)]
    }
  # precompute quantities used later
  logfactorialy <- lgamma(y+1)
  log_lambda <- log(lambda)
  log.Z <- logZ(log_lambda, nu, summax = summax)
  # compute loglikelihood 
  loglik <- sum(y*log_lambda - nu*logfactorialy - log.Z)
  return(loglik)
}

#' @rdname comp_mu_loglik
comp_mu_neg_loglik_log_nu_only <-function(log_nu, mu, y, summax){
  # compute negative loglikelihood for COMP-mu regression models
  # precompute quantities used later
  nu = exp(log_nu)
  logfactorialy <- lgamma(y+1)
  lambda <- (mu+(nu-1)/(2*nu))^(nu)
  log_lambda <- log(lambda)
  log.Z <- logZ(log_lambda, nu, summax)
  #meanlogfactorialy <- comp_mean_logfactorialy(lambda, nu, mu)
  # compute loglikelihood 
  nloglik <- -sum(y*log_lambda - nu*logfactorialy - log.Z)
  return(nloglik)
}
