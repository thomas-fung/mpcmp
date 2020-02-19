#' Functions to Compute Various Expected Values for the COM-Poisson Distribution
#' 
#' Functions to approximate the various expected values for the COM-Poisson distribution
#' via truncation.  The standard COM-Poisson parametrization is being used here. 
#' The lambda and nu values are recycled to match the length 
#' of the longer one and that would determine the length of the results. 
#' @param lambda,nu, rate and dispersion  parameters. Must be positives. 
#' @param log.Z, an optional vector specifying normalizing constant Z in log scale. 
#' @param summax maximum number of terms to be considered in the truncated sum. The
#' default is to sum to 100. 
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
#' @name comp_expected_values
NULL

#' @rdname comp_expected_values
comp_mean_logfactorialy = function(lambda, nu, log.Z, summax=100){
  # approximates mean by truncation of Ylog(Y!) for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  if (missing(log.Z)){
    df <- CBIND(lambda=lambda, nu=nu)
    log.Z <- logZ(log(df[,1]), df[,2], summax)
  }
  df <- CBIND(lambda=lambda, nu=nu, log.Z)
  lambda <- df[,1]
  nu <- df[,2]
  log.Z <- df[,3]
  term <- matrix(0, nrow = length(lambda), ncol=summax)
  for (y in 1:summax){
    term[,y] <- lgamma(y)*exp((y-1)*log(lambda) - nu*lgamma(y)-log.Z)
  }
  mean1 <- apply(term,1,sum)
  return(mean1)
}

#' @rdname comp_expected_values
#' @export
comp_mean_ylogfactorialy <- function(lambda, nu, log.Z, summax=100){
  # approximates mean by truncation of Ylog(Y!) for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  if (missing(log.Z)){
    df <- CBIND(lambda=lambda, nu=nu)
    log.Z <- logZ(log(df[,1]), df[,2], summax)
  }
  df <- CBIND(lambda=lambda, nu=nu, log.Z = log.Z)
  lambda <- df[,1]
  nu <- df[,2]
  log.Z <- df[,3]
  term <- matrix(0, nrow = length(lambda), ncol=summax)
  for (y in 1:summax) {
    term[,y] <- exp(log(y-1)+log(lgamma(y))+ (y-1)*log(lambda) - nu*lgamma(y)-log.Z)
  }
  mean1 <- apply(term,1,sum)
  return(mean1)
}

#' @rdname comp_expected_values
#' @export
comp_means <- function(lambda, nu, log.Z, summax=100) {
  # approximates mean by truncation of COMP distributions
  # lambda, nu, mu.bd are recycled to match the length of each other.
  if (missing(log.Z)){
    df <- CBIND(lambda=lambda, nu=nu)
    log.Z <- logZ(log(df[,1]), df[,2], summax)
  }
  df <- CBIND(lambda=lambda, nu=nu, log.Z = log.Z)
  lambda <- df[,1]
  nu <- df[,2]
  log.Z <- df[,3]
  term <- matrix(0, nrow = length(lambda), ncol=summax)
  for (y in 1:summax) {
    term[,y] <- exp(log(y-1)+(y-1)*log(lambda) - nu*lgamma(y)-log.Z)
  }
  mean1 <- apply(term,1,sum)
  return(mean1)
}


#' @rdname comp_expected_values
#' @export
comp_variances <- function(lambda, nu, log.Z, summax=100) {
  # approximates normalizing constant by truncation for COMP distributions
  # lambda, nu, mu.bd are recycled to match the length of each other.
  if (missing(log.Z)){
    df <- CBIND(lambda=lambda, nu=nu)
    log.Z <- logZ(log(df[,1]), df[,2], summax)
  }
  df <- CBIND(lambda=lambda, nu=nu, log.Z = log.Z)
  lambda <- df[,1]
  nu <- df[,2]
  log.Z <- df[,3]
  term <- matrix(0, nrow = length(lambda), ncol=summax)
  for (y in 1:summax){
    term[,y] <- exp(2*log(y-1) + (y-1)*log(lambda) - nu*lgamma(y)-log.Z)
  }
  var1 <- apply(term,1,sum)- (comp_means(lambda, nu, log.Z, summax))^2
  return(var1)
}

#' @rdname comp_expected_values
#' @export
comp_variances_logfactorialy <- function(lambda, nu, log.Z, summax = 100) {
  # approximates normalizing constant by truncation for COMP distributions
  # lambda, nu, mu.bd are recycled to match the length of each other.
  if (missing(log.Z)){
    df <- CBIND(lambda=lambda, nu=nu)
    log.Z <- logZ(log(df[,1]), df[,2], summax)
  }
  df <- CBIND(lambda=lambda, nu=nu, log.Z = log.Z)
  lambda <- df[,1]
  nu <- df[,2]
  log.Z <- df[,3]
  term <- matrix(0, nrow = length(lambda), ncol=summax)
  for (y in 1:summax){
    term[,y] <- lgamma(y)^2*exp((y-1)*log(lambda) - nu*lgamma(y)- log.Z)
  }
  var1 <- apply(term,1,sum) - (comp_mean_logfactorialy(lambda, nu, log.Z, summax))^2
  return(var1)
}
