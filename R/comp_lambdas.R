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
#' @param summax maximum number of terms to be considered in the truncated sum
#' @return 
#' Both \code{comp_lambdas} and \code{comp_lambdas_fixed_ub} returns the lambda value(s)
#' that satisfies the mean constraint(s) as well as the current lambdaub value. 
#' lambda value(s) returns by \code{comp_lambdas_fixed_ub} is bounded by the lambdaub 
#' value. 
#' \code{comp_lambdas} has the extra ability to scale up/down lambdaub to find the most 
#' appropriate lambda values. 
#'  
#' @name comp_lambdas
NULL

#' @rdname comp_lambdas
comp_lambdas <- function(mu, nu, lambdalb = 1e-10, lambdaub = 1000, 
                         maxlambdaiter = 1e3, tol = 1e-6, lambdaint = 1, summax = 100){
  df <- CBIND(mu=mu, nu=nu, lambda = lambdaint, lb = lambdalb, ub = lambdaub)
  mu <- df[,1]
  nu <- df[,2]
  lambda <- df[,3]
  lambdalb <- df[,4]
  lambdaub <- df[,5]
  lambda.ok <- 
    comp_lambdas_fixed_ub(mu, nu, lambdalb = lambdalb, lambdaub = lambdaub, 
                          maxlambdaiter = maxlambdaiter, tol = tol, 
                          lambdaint = lambda, summax = summax)
  lambda <- lambda.ok$lambda
  lambdaub <- lambda.ok$lambdaub
  lambdaub.err.ind <- which(is.nan(lambda))
  sub_iter1 <- 1
  while (length(lambdaub.err.ind)>0 && sub_iter1 <=100){
    lambdaub[lambdaub.err.ind] <- 0.5*lambdaub[lambdaub.err.ind]
    lambda.ok <- 
      comp_lambdas_fixed_ub(mu[lambdaub.err.ind], nu[lambdaub.err.ind], 
                            lambdalb = lambdalb[lambdaub.err.ind], 
                            lambdaub = lambdaub[lambdaub.err.ind], 
                            maxlambdaiter = maxlambdaiter, tol = tol, 
                            lambdaint = lambda[lambdaub.err.ind], summax = summax)
    lambda[lambdaub.err.ind] <- lambda.ok$lambda
    lambdaub[lambdaub.err.ind] <- lambda.ok$lambdaub
    sub_iter1 <- sub_iter1+1
    lambdaub.err.ind <-  which(is.nan(lambda))
  }
  lambdaub.err.ind <- which(lambda/lambdaub >= 1-tol)
  sub_iter1 <- 1 
  while (length(lambdaub.err.ind)>0 && sub_iter1 <= 100){
    lambdaub[lambdaub.err.ind] <- 2^(sub_iter1)*lambdaub[lambdaub.err.ind]
    lambda.ok <- 
      comp_lambdas_fixed_ub(mu[lambdaub.err.ind], nu[lambdaub.err.ind], 
                            lambdalb = lambdalb[lambdaub.err.ind], 
                            lambdaub = lambdaub[lambdaub.err.ind], 
                            maxlambdaiter = maxlambdaiter, tol = tol, 
                            lambdaint = lambda[lambdaub.err.ind],
                            summax = summax)
    lambda[lambdaub.err.ind] <- lambda.ok$lambda
    lambdaub[lambdaub.err.ind] <- lambda.ok$lambdaub
    sub_iter1 <- sub_iter1+1
    lambdaub.err.ind <- which(lambda/lambdaub >= 1-tol)
  }
  out <- list()
  out$lambda <- lambda
  out$lambdaub <- max(lambdaub)
  return(out)
}

#' @rdname comp_lambdas
comp_lambdas_fixed_ub <- function(mu, nu, lambdalb = 1e-10, lambdaub = 1000, 
                                  maxlambdaiter = 1e3, tol = 1e-6, lambdaint = 1, 
                                  summax = 100) {
  #df <- CBIND(mu=mu, nu=nu, lambda = lambdaint, lb = lambdalb, ub = lambdaub)
  #mu <- df[,1]
  #nu <- df[,2]
  #lambda <- df[,3]
  #lb <- df[,4]
  #ub <- df[,5]
  lambda <- lambdaint
  lb <- lambdalb
  ub <- lambdaub
  iter <- 1
  log.Z <- logZ(log(lambda), nu, summax = summax)
  mean1 <- comp_means(lambda, nu, log.Z = log.Z, summax = summax)
  not.converge.ind <- which(abs(mean1-mu)>tol)
  while (length(not.converge.ind)>0 && iter <200){
    still.above.target.ind <-  which((mean1[not.converge.ind]
                                    >mu[not.converge.ind]))
    still.below.target.ind <-  which(mean1[not.converge.ind] < mu[not.converge.ind])
    lb[not.converge.ind[still.below.target.ind]] <- 
      lambda[not.converge.ind[still.below.target.ind]]
    ub[not.converge.ind[still.above.target.ind]] <- 
      lambda[not.converge.ind[still.above.target.ind]]
    lambda <-  (lb+ub)/2
    log.Z <- logZ(log(lambda), nu, summax = summax)
    mean1 <- comp_means(lambda, nu, log.Z = log.Z, summax = summax)
    while (sum(mean1==0)>0){
      ub[not.converge.ind[mean1==0]] <- ub[not.converge.ind[mean1==0]]/2
      lambdaub <-  lambdaub/2
      lambda <-  (lb+ub)/2
      log.Z <- logZ(log(lambda), nu, summax = summax)
      mean1 <- comp_means(lambda, nu, log.Z = log.Z, summax = summax)
    }
    not.converge.ind <- which((1-(((abs(mean1-mu) <=tol) + (lambda == lb) + (lambda == ub)
                                   + (ub == lb)) >= 1))==1)
    iter <- iter+1
  }
  while (length(not.converge.ind)>0 && iter <maxlambdaiter){
    # basically the content of comp_variances without recalculating Z and mean1
    term <- matrix(0, nrow = length(mu), ncol=summax)
    for (y in 1:summax){
      term[,y] <- exp(2*log(y-1)+(y-1)*log(lambda) - nu*lgamma(y)-log.Z)
    }
    var1 <- apply(term,1,sum) - mean1^2
    ## newton raphson update
    newtonsteps <- - lambda[not.converge.ind]*mean1[not.converge.ind]/
      (var1[not.converge.ind])^2*(log(mean1[not.converge.ind])-log(mu[not.converge.ind]))
    
    
    lambda.new <- lambda[not.converge.ind] + newtonsteps
    ## if newton raphson steps out of bound, use bisection method
    out.of.bound.ind <-  which((lambda.new< lb[not.converge.ind])
                             + (lambda.new>ub[not.converge.ind]) ==1)
    if (length(out.of.bound.ind>0)){
      lambda.new[out.of.bound.ind] <- 
        (lb[not.converge.ind[out.of.bound.ind]]+ub[not.converge.ind[out.of.bound.ind]])/2
      # any out of bound updates are replaced with mid-point of ub and lb
    }
    lambda[not.converge.ind] <- lambda.new
    log.Z <- logZ(log(lambda), nu, summax)
    term <- matrix(0, nrow = length(mu), ncol=summax)
    for (y in 1:summax) {
      term[,y] <- exp(log(y-1)+(y-1)*log(lambda) - nu*lgamma(y)-log.Z)
    }
    mean1 <- apply(term,1,sum)
    if (length(out.of.bound.ind)>0){
      still.above.target.ind <-  which(mean1[not.converge.ind[out.of.bound.ind]]
                                     >mu[not.converge.ind[out.of.bound.ind]])
      still.below.target.ind <-  which(mean1[out.of.bound.ind]<mu[out.of.bound.ind])
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
    not.converge.ind <- which((1-(((abs(mean1-mu) <=tol) + (lambda == lb) + (lambda == ub)
                                   + (ub ==lb)) >= 1))==1)
    iter <- iter+1
  }
  out <- list()
  out$lambda <- lambda
  out$lambdaub <- lambdaub
  return(out)
}


