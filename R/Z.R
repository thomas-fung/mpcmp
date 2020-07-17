#' Calculate the Normalizing Constant for COM-Poisson distribution
#' 
#' A function to approximate the normalizing constant for COM-Poisson distributions
#' via truncation. The standard COM-Poisson parametrization is being used here. 
#' Notice that the sum is hard coded to truncate at 100 so the approximation will be quite
#' bad if the COM-Poisson has a large rate or mean. 
#' 
#' @param lambda rate parameter, straightly positive
#' @param nu dispersion parameter, straightly positive
#' @param log.z logical; if \code{TRUE}, normalising constant \eqn{Z} are returned as 
#' \eqn{log(Z)}.
#' @param summax maximum number of terms to be considered in the truncated sum
Z <- function(lambda, nu, log.z = FALSE, summax){
  # approximates normalizing constant for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  .Deprecated("logZ")
  df <- CBIND(lambda=lambda, nu=nu)
  lambda <- df[,1]
  nu <- df[,2]
  termlim <- 1e-6
  # zero vector length of same length as lambda
  sum1 <- term <- rep(0, length(lambda))
  not.converge.ind <- 1:length(lambda)
  for(y in 1:summax) {
    term[not.converge.ind] <- exp((y-1)*log(lambda[not.converge.ind]) - 
                                      nu[not.converge.ind]*lgamma(y))
    sum1[is.infinite(term)] <- Inf
    term[is.infinite(sum1)] <- 0
    if (y > 3){
      not.converge.ind = which(term/sum1 >= termlim)
      if (length(not.converge.ind)==0) {
        break
      }
    }
    sum1 <- sum1+term
  }
  sum1 <- log(sum1)
  infinite.ind <- which(is.infinite(sum1))
  if (length(infinite.ind)!=0){
    nu.sub <- nu[infinite.ind]
    lambda.sub <- lambda[infinite.ind]
    A <- (8*nu.sub^2+12*nu.sub+3)/(96*(lambda.sub^(1/nu.sub))*nu.sub^2)
    B <- (1+6*nu.sub)/(144*(lambda.sub^(2/nu.sub))*nu.sub^3)
    sum1.sub <- nu.sub*lambda.sub^(1/nu.sub)-
      ((nu.sub-1)/(2*nu.sub)*log(lambda.sub)+(nu.sub-1)/2*log(2*pi)+0.5*log(nu.sub))+
      log(1+(nu.sub-1)*(A+B))
      sum1[infinite.ind] <- sum1.sub
  }
  if (!log.z){
    sum1 <- exp(sum1)}
  return(sum1)
}



#' Calculate the Normalizing Constant in log scale for COM-Poisson distribution
#' 
#' A function to approximate the normalizing constant for COM-Poisson distributions via 
#' truncation. The standard COM-Poisson parametrization is being used here. 
#' 
#' As of version 0.2.0 of this package, \code{logZ} will supersede \code{Z} for calculating
#' the normalizing constant. \code{logZ} utilised a method that can calculate
#'  \code{log(exp(logx) + exp(logy))} in a somewhat numerically stable way. 
#' 
#' This function was originally purposed in the \code{cmpreg} package of Ribeiro Jr, 
#' Zeviani & Demétrio (2019).
#' 
#' @param log_lambda rate parameter in log scale.
#' @param nu dispersion parameter, straightly positive.
#' @param summax maximum number of terms to be considered in the truncated sum.
#' @references 
#'  Ribeiro Jr, E. E., Zeviani, W. M., Demétrio, C. G. B. (2019) \code{cmpreg}: 
#'  Reparametrized COM-Poisson Regression Models. R package version 0.0.1.

logZ <- function(log_lambda, nu, summax = 100){
  # approximates normalizing constant for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(log_lambda=log_lambda, nu=nu)
  log_lambda <- df[,1]
  nu <- df[,2]
  return(logZ_c(log_lambda,nu,summax=summax))
  }

