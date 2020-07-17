#' Parameter Generator for nu
#' 
#' A function that use the arguments of a \code{glm.cmp} call to generate a better initial 
#' \code{nu} estimate. 
#' 
#' From version 0.3.4, this function is no longer being used as part of the estimation algorithm 
#' and this function will be defunct in our next update.
#' 
#' @param param numeric vector: the model coefficients & the current value of \code{nu}. 
#' It is assumed that \code{nu} is in the last position of \code{param}.
#' @param y numeric vector: response variable
#' @param xx numeric matrix: the explanatory variables
#' @param offset numeric vector: a vector of length equal to the number of cases
#' @param llstart numeric: current log-likelihood value
#' @param fsscale numeric: a scaling factor (generally >1) for the 
#' relaxed fisher scoring algorithm
#' @param lambdalb,lambdaub numeric: the lower and upper end points for the interval to be
#' searched for lambda(s). 
#' @param maxlambdaiter numeric: the maximum number of iterations allowed to solve 
#' for lambda(s).
#' @param tol numeric: the convergence threshold. A lambda is said to satisfy the 
#' mean constraint if the absolute difference between the calculated mean and a fitted
#' values is less than tol.
#' @param summax maximum number of terms to be considered in the truncated sum
#' @return 
#' List containing the following:
#' \item{param}{the model coefficients & the updated \code{nu}}
#' \item{maxl}{the updated log-likelihood}
#' \item{fsscale}{the final scaling factor used}
#' 
getnu <- function(param, y, xx , offset, llstart, fsscale = 1, 
                  lambdalb = 1e-10, lambdaub = 1000, maxlambdaiter = 1e3, tol = 1e-6,
                  summax  = 100){
  options(warn=2)
  n <- length(y) # sample size
  q <- ncol(xx)  # number of covariates
  nu_lb <- 1e-10
  iter <- 1
  beta <- param[1:q]
  mu <-  exp(t(xx%*%beta)[1,]+offset)
  lambda <- param[(q+1):(q+n)]
  nu_old <- param[q+n+1]
  ll_old <- llstart
  log.Z <- logZ(log(lambda), nu_old, summax = summax)
  ylogfactorialy <- comp_mean_ylogfactorialy(lambda, nu_old, log.Z, summax)
  logfactorialy <- comp_mean_logfactorialy(lambda, nu_old, log.Z, summax)
  variances <- comp_variances(lambda, nu_old, log.Z, summax)
  variances_logfactorialy <- comp_variances_logfactorialy(lambda,nu_old, log.Z, summax)
  Aterm <- (ylogfactorialy- mu*logfactorialy)
  update_score <- sum(Aterm*(y-mu)/variances -(lgamma(y+1)-logfactorialy))
  update_info_matrix <- sum(-Aterm^2/variances+ variances_logfactorialy)
  nu <- nu_old + fsscale*update_score/update_info_matrix
  while (nu < nu_lb){
    nu <- (nu+nu_old)/2
  }
  lambdaold <- lambda
  lambda.ok <- comp_lambdas(mu, nu, lambdalb = lambdalb, 
                            #lambdaub = lambdaub,
                            lambdaub = min(lambdaub,2*max(lambdaold)),
                            maxlambdaiter = maxlambdaiter, tol = tol, 
                            lambdaint = lambdaold, summax = summax)
  lambda <- lambda.ok$lambda
  lambdaub <- lambda.ok$lambdaub
  param <- c(beta, lambda, nu)
  ll_new <- comp_mu_loglik(param = param, y=y, xx= xx, offset= offset, summax= summax)
  sub_iter <- 1
  while (ll_new < ll_old && sub_iter < 20 && abs((ll_new-ll_old)/ll_new)>tol){
    nu <- (nu+nu_old)/2
    fsscale <- max(fsscale/2,1)
    lambda.ok <- comp_lambdas(mu, nu, lambdalb = lambdalb, 
                              #lambdaub = lambdaub,
                              lambdaub = min(lambdaub,max(2*lambdaold)), 
                              maxlambdaiter = maxlambdaiter, tol = tol, 
                              lambdaint=lambda, summax =summax)
    lambda <- lambda.ok$lambda
    lambdaub<- lambda.ok$lambdaub
    param <- c(beta, lambda, nu)
    ll_new <- comp_mu_loglik(param = param, y=y, xx= xx, offset= offset, summax= summax)
    sub_iter <- sub_iter + 1
  }
  iter = iter+1
  while ((abs((ll_new-ll_old)/ll_new)>1e-6) && abs(nu-nu_old)>1e-6 && iter <=100){
    ll_old <- ll_new
    nu_old <- nu
    log.Z = logZ(log(lambda),nu_old, summax = summax)
    ylogfactorialy <- comp_mean_ylogfactorialy(lambda,nu_old, log.Z, summax)
    logfactorialy <- comp_mean_logfactorialy(lambda,nu_old, log.Z, summax)
    variances <- comp_variances(lambda,nu_old, log.Z, summax)
    variances_logfactorialy <- comp_variances_logfactorialy(lambda,nu_old, log.Z, summax)
    Aterm <- (ylogfactorialy- mu*logfactorialy)
    update_score <- sum(Aterm*(y-mu)/variances -(lgamma(y+1)-logfactorialy))
    update_info_matrix <- sum(-Aterm^2/variances+ variances_logfactorialy)
    nu <- nu_old + fsscale*update_score/update_info_matrix
    while (nu < nu_lb){
      nu <- (nu + nu_old)/2
    }
    if (abs(nu-nu_old)<1e-6){
      break
    }
    lambdaold <- lambda
    lambda.ok <- comp_lambdas(mu,nu, lambdalb = lambdalb, 
                              #lambdaub = lambdaub,
                              lambdaub = min(lambdaub,2*max(lambdaold)), 
                              maxlambdaiter = maxlambdaiter, tol = tol, 
                              lambdaint= lambda, summax= summax)
    lambda <- lambda.ok$lambda
    lambdaub <- lambda.ok$lambdaub
    param <- c(beta, lambda, nu)
    ll_new <- comp_mu_loglik(param = param, y=y, xx= xx, offset= offset, summax = summax)
    sub_iter <- 1
    while (ll_new < ll_old && sub_iter < 20 && abs((ll_new-ll_old)/ll_new)>tol){
      nu <- (nu+nu_old)/2
      fsscale <- max(fsscale/2,1)
      lambdaold <- lambda
      lambda.ok <- comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                #lambdaub = lambdaub, 
                                lambdaub = min(lambdaub,2*max(lambdaold)), 
                                maxlambdaiter = maxlambdaiter, tol = tol, 
                                lambdaint = lambda, summax = summax)
      lambda <- lambda.ok$lambda
      lambdaub <- lambda.ok$lambdaub
      param <- c(beta, lambda, nu)
      ll_new <- comp_mu_loglik(param = param, y=y, xx= xx, offset= offset, summax=summax)
      sub_iter <- sub_iter+1
    }
    iter = iter+1
  }
  obj <- list()
  obj$param <- param
  obj$maxl <- ll_new
  obj$fsscale <- fsscale
  obj$iter <- iter
  obj$lambdaub <- lambdaub
  options(warn=0)
  return(obj)
}
