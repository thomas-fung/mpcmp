#' Functions to Compute Various Expected Values for the COM-Poisson Distribution
#' 
#' Functions to approximate the various expected values for the COM-Poisson distribution
#' via truncation.  The standard COM-Poisson parametrization is being used here. 
#' The lambda and nu values are recycled to match the length 
#' of the longer one and that would determine the length of the results. 
#' Notice that the sum is hard coded to tuncate at 100 so the approximation will be quite
#' bad if the COM-Poisson has a large rate or mean. 
#' 
#' @param lambda,nu rate and dispersion parameters. Must be positives. 
#' @return 

#' \code{comp_mean_logfactorialy} gives the mean of \emph{log(Y!)}. 
#' 
#' \code{comp_mean_ylogfactorialy} gives the mean of \emph{ylog(Y!)}.
#' 
#' \code{comp_means} gives the mean of \emph{Y}. 
#' 
#' \code{comp_variances} gives the variance of \emph{Y}.
#' 
#' \code{comp_variances_logfactorialy} gives the variance of \emph{log(Y!)}.
#' 
#' @name COMP_Expected_Values
NULL

#' @rdname COMP_Expected_Values
comp_mean_logfactorialy = function(lambda, nu){
  # approximates mean by truncation of Ylog(Y!) for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(lambda=lambda, nu=nu)
  lambda <- df[,1]
  nu <- df[,2]
  if (length(lambda)>1 && length(nu>1) && length(lambda)!= length(nu)){
    stop("lambda, nu must be scalars or vectors of the same length")}
  summax <- 100
  termlim <- 1e-6
  sum1 <- 0
  for (y in 1:summax){
    term <- lgamma(y)*exp(log(lambda^(y-1)) - nu*lgamma(y))
    if (y > 3) {
      if (max(term/sum1, na.rm = TRUE) < termlim){
        break
      }
    }
    sum1 <- sum1 + term
  }
  mean1 <- sum1/Z(lambda, nu)
  return(mean1)
}

#' @rdname COMP_Expected_Values
comp_mean_ylogfactorialy <- function(lambda, nu){
  # approximates mean by truncation of Ylog(Y!) for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(lambda=lambda, nu=nu)
  lambda <- df[,1]
  nu <- df[,2]
  summax <- 100
  termlim <- 1e-6
  sum1 <- 0
  for (y in 1:summax) {
    term <- (y-1)*lgamma(y)*exp(log(lambda^(y-1)) - nu*lgamma(y))
    if (y > 3) {
      if ( max(term/sum1,na.rm = TRUE) < termlim){
        break
      }
    }
    sum1 <- sum1 + term
  }
  mean1 <- sum1/Z(lambda, nu)
  return(mean1)
}

#' @rdname COMP_Expected_Values
comp_means <- function(lambda, nu) {
  # approximates mean by truncation of COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(lambda=lambda, nu=nu)
  lambda <- df[,1]
  nu <- df[,2]
  summax <- 100
  termlim <- 1e-6
  sum1 <- 0
  for (y in 1:summax) {
    term <- (y-1)*exp(log(lambda^(y-1)) - nu*lgamma(y))
    if (y > 3) {
      if (max(term/sum1,na.rm = TRUE) < termlim){
        break
      }
    }
    sum1 <- sum1 + term
    mean1 <- sum1/Z(lambda, nu)
  }
  return(mean1)
}

#' @rdname COMP_Expected_Values
comp_variances <- function(lambda, nu) {
  # approximates normalizing constant by truncation for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(lambda=lambda, nu=nu)
  lambda <- df[,1]
  nu <- df[,2]
  summax <- 100
  termlim <- 1e-6
  sum2 <- 0
  for (y in 1:summax){
    term <- (y-1)^2*exp(log(lambda^(y-1)) - nu*lgamma(y))
    if (y > 3) {
      if (max(term/sum2,na.rm = TRUE) < termlim) {
        break
      }
    }
    sum2 <- sum2 + term
  }
  var1 <- sum2/Z(lambda,nu) - (comp_means(lambda,nu))^2
  return(var1)
}

#' @rdname COMP_Expected_Values
comp_variances_logfactorialy <- function(lambda, nu) {
  # approximates normalizing constant by truncation for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(lambda=lambda, nu=nu)
  lambda <- df[,1]
  nu <- df[,2]
  summax <- 100
  termlim <- 1e-6
  sum2 <- 0
  for (y in 1:summax){
    term <- (lgamma(y))^2*exp(log(lambda^(y-1)) - nu*lgamma(y))
    if (y > 3) {
      if (max(term/sum2,na.rm = TRUE) < termlim) {
        break
      }
    }
    sum2 <- sum2 + term
  }
  var1 <- sum2/Z(lambda,nu) - (comp_mean_logfactorialy(lambda,nu))^2
  return(var1)
}
