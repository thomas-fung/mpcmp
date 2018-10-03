#' Solve for Lambda for a Particular Mean Parametrized COM-Poisson Distribution
#' 
#' Given a particular mean parametrized COM-Poisson distribution i.e. mu and nu, 
#' this function is used to find a lambda that can satisfy the mean constraint with a 
#' combination of bisection and Newton-Raphson updates. The function is also vectorized but
#' will only update those that have not converged. 
#' 
#' @param mu,nu mean and dispersion parameters. Must be straightly positive. 
#' @param lambdalb,lambdaub numeric; the lower and upper end points for the interval to be
#' searched for lambda(s). 
#' @param maxlambdaiter numeric; the maximum number of iterations allowed to solve 
#' for lambda(s).
#' @param tol numeric; the convergence threshold. A lambda is said to satisfy the 
#' mean constraint if the absolute difference between the calculated mean and the 
#' corresponding mu values is less than tol. 
#' @param lambdaint numeric vector; initial gauss for lambda(s).
#' @return 
#' The function returns the lambda value(s) that satisfies the mean constraint(s). 
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