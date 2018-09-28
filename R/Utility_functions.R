#' Combine R Objects by Columns
#' 
#' Take a sequence of vector, matrix or data-frame arguments and combine them by columns. 
#' \code{CBIND} is used within the package over \code{cbind} to recycle the shorter 
#' arguments so that their number of rows would match. 
#' 
#' @param ... (generalized) vectors or matrices. These can be given as named arugments. 
#' @param deparse.level integer; deparse.level = 0 constructs no labels, 
#' deparse.level = 1 (the default) or > 1 constructs labels from the arguments names. 
#' 
CBIND <- function(..., deparse.level = 1) {
  dots <- list(...)
  len <- sapply(dots, length)
  dots <- lapply(seq_along(dots),
                 function(i, x, len) rep(x[[i]], length.out = len),
                 x = dots, len = max(len))
  do.call(cbind, c(dots, deparse.level = deparse.level))
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}


Z <- function(lambda, nu){
  # approximates normalizing constant for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(lambda=lambda, nu=nu)
  lambda <- df[,1]
  nu <- df[,2]
  summax <- 100
  termlim <- 1e-6
  # zero vector length of same length as lambda
  sum1 <- 0
  for(y in 1:summax) {
    term <- exp(log(lambda^(y-1)) - nu*lgamma(y))
    if (y > 3){
      if (max(term/sum1,na.rm = TRUE) < termlim) {
        break
      }
    }
    sum1 <- sum1 + term
  }
  return(sum1)
}


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


comp_lambdas <- function(mu, nu, lambdalb = 1e-10, lambdaub = 1900,
                         maxlambdaiter = 1e3, tol = 1e-6, lambdaint = 1) {
  # for given vector mean mu and vector dispersion nu, solve for the vector rate lambda
  df <- CBIND(mu=mu, nu=nu, lambda = lambdaint, lb = lambdalb, ub = lambdaub)
  mu <- df[,1]
  nu <- df[,2]
  lambda <- df[,3]
  lb <- df[,4]
  ub <- df[,5]
  # basically the contecnt of comp_means so that we don't have to
  # recalculate Z within an iteration
  summax <- 100
  termlim <- 1e-6
  iter <- 1
  Zcurrent <- Z(lambda, nu)
  sum1 <- 0
  for (y in 1:summax) {
    term <- (y-1)*exp(log(lambda^(y-1)) - nu*lgamma(y))
    if (y > 3) {
      if (max(term/sum1,na.rm = TRUE) < termlim){
        break
      }
    }
    sum1 <- sum1 + term
  }
  mean1 <- sum1/Zcurrent
  ### end of comp_means function
  not.converge.ind = which(abs(mean1-mu)>tol)
  while (length(not.converge.ind)>0 && iter <200){
    still.above.target.ind = which(mean1[not.converge.ind]
                                   >mu[not.converge.ind])
    still.below.target.ind = which(mean1[not.converge.ind]<mu[not.converge.ind])
    lb[not.converge.ind[still.below.target.ind]] =
      lambda[not.converge.ind[still.below.target.ind]]
    ub[not.converge.ind[still.above.target.ind]] =
      lambda[not.converge.ind[still.above.target.ind]]
    lambda = (lb+ub)/2
    Zcurrent <- Z(lambda, nu)
    sum1 <- 0
    for (y in 1:summax) {
      term <- (y-1)*exp(log(lambda^(y-1)) - nu*lgamma(y))
      if (y > 3) {
        if (max(term/sum1,na.rm = TRUE) < termlim){
          break
        }
      }
      sum1 <- sum1 + term
    }
    mean1 <- sum1/Zcurrent
    not.converge.ind <- which(abs(mean1-mu)>tol)
    iter <- iter+1
  }
  while (length(not.converge.ind)>0 && iter <maxlambdaiter){
    # basically the content of comp_variances without recalculating Z
    sum2 <- 0;
    for (y in 1:summax){
      term <- (y-1)^2*exp(log(lambda^(y-1)) - nu*lgamma(y))
      if (y > 3) {
        if (max(term/sum2,na.rm = TRUE) < termlim) {
          break
        }
      }
      sum2 <- sum2 + term
    }
    var1 <- sum2/Zcurrent - mean1^2
    ## newton raphson update
    newtonsteps <- - lambda[not.converge.ind]*mean1[not.converge.ind]/(var1[not.converge.ind])^2*
      (log(mean1[not.converge.ind])-log(mu[not.converge.ind]))
    lambda.new <- lambda[not.converge.ind] + newtonsteps
    ## if newton raphson steps out of bound, use bisection method
    out.of.bound.ind = which((lambda.new< lb[not.converge.ind])
                             + (lambda.new>ub[not.converge.ind]) ==1)
    if (length(out.of.bound.ind>0)){
      lambda.new[out.of.bound.ind] =
        (lb[not.converge.ind[out.of.bound.ind]]+ub[not.converge.ind[out.of.bound.ind]])/2
      # any out of bound updates are replaced with mid-point of ub and lb
    }
    lambda[not.converge.ind] <- lambda.new
    Zcurrent <- Z(lambda, nu)
    sum1 <- 0
    for (y in 1:summax) {
      term <- (y-1)*exp(log(lambda^(y-1)) - nu*lgamma(y))
      if (y > 3) {
        if (max(term/sum1,na.rm = TRUE) < termlim){
          break
        }
      }
      sum1 <- sum1 + term
    }
    mean1 <- sum1/Zcurrent
    if (length(out.of.bound.ind)>0){
      still.above.target.ind = which(mean1[not.converge.ind[out.of.bound.ind]]
                                     >mu[not.converge.ind[out.of.bound.ind]])
      still.below.target.ind = which(mean1[out.of.bound.ind]<mu[out.of.bound.ind])
      if (length(still.below.target.ind)>0){
        lb[not.converge.ind[out.of.bound.ind[still.below.target.ind]]] =
          lambda[not.converge.ind[out.of.bound.ind[still.below.target.ind]]]
      }
      if (length(still.above.target.ind)>0){
        ub[not.converge.ind[out.of.bound.ind[still.above.target.ind]]] =
          lambda[not.converge.ind[out.of.bound.ind[still.above.target.ind]]]
        # any out of bound updates are replaced with mid-point of ub and lb
      }
    }
    not.converge.ind <- which(((abs(mean1-mu)>tol)+ (lambda != lb)+ (lambda!=ub)) >= 1)
    iter <- iter+1
  }
  return(lambda)
}
