
#' Parameter Generator for nu
#' 
#' A function that use the arguments of a \code{glm.cmp.gamma} call to generate a better initial 
#' \code{nu} estimate. 
#' 
#' @param param numeric vector: the model coefficients & the current value of \code{nu}. 
#' @param y numeric vector: response variable
#' @param xx numeric matrix: the explanatory variables
#' @param offset numeric vector: a vector of lenght equal to the number of cases
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

getnu.gamma<- function(param, y, xx=X , S = S, offset, llstart, fsscale = 1, lambdalb = 1e-10, 
                       lambdaub = 1900, maxlambdaiter = 1000, tol = 1e-06) {
  options(warn=2)
  summax  <- 300
  n <- length(y) 
  q <- ncol(xx) 
  t <- ncol(S)                       
  nu_lb <- 1e-10
  iter <- 1
  beta <- param[1:q]
  mu <-  exp(t(xx%*%beta)[1,]) 
  gam_old <- param[(q+2*n+1):(q+2*n+t)]
  nu_old <- param[(q+n+1):(q+2*n)]        
  ll_old <- llstart  #opposit
  lambda <- param[(q+1):(q+n)]                   
  update_score <- matrix(0,ncol=t)
  update_info_matrix <-matrix(0,ncol=t,nrow=t)
  log.Z <- logZ(log(lambda), nu_old, summax = summax)
  Aterm <- (comp_mean_ylogfactorialy(lambda, nu_old, log.Z, summax = summax)-                  
              mu*comp_mean_logfactorialy(lambda,nu_old, log.Z, summax = summax)) 
  for (k in 1:n){
    update_score <- update_score+ sum(Aterm[k] * (y[k] - mu[k])/comp_variances(lambda[k], nu_old[k], log.Z[k], summax = summax)
                                      - (lgamma(y[k] + 1) - comp_mean_logfactorialy(lambda[k], nu_old[k], log.Z[k], summax = summax)))*(nu_old[k]%*%S[k,]) 
    update_info_matrix <- update_info_matrix + sum(-Aterm[k]^2/comp_variances(lambda[k], nu_old[k], log.Z[k], summax = summax)
                                                   +  comp_variances_logfactorialy(lambda[k], nu_old[k], log.Z[k], summax = summax))*(nu_old[k]^2*S[k,]%*%t(S[k,]))}
  
  gam <- gam_old + fsscale*update_score%*%solve(update_info_matrix)
  eta1 <- t(S%*%t(gam))[1,]
  nu <- exp(eta1)
  
  
  lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, lambdaub = lambdaub, 
                             maxlambdaiter = maxlambdaiter, tol = tol, summax = summax))
  
  
  if (class(lambda) == "try-error") {
    while (class(lambda) == "try-error") {
      lambdaold <- NULL
      lambdaubold <- lambdaub
      lambdaub <- 0.8 * lambdaub
      lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                 lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                 tol = tol))
    }
    lastworking_lambdaub <- lambdaub
    sub_iter1 <- 1
    while (max(lambda$lambda)/max(lambdaub) >= 1 - tol && sub_iter1 <= 
           20) {
      lambdaubnew <- lambdaub
      lambdaub <- (lambdaubnew + lambdaubold)/2
      sub_iter1 <- sub_iter1 + 1
      lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                 lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                 tol = tol))
      sub_iter2 <- 1
      while (class(lambda) == "try-error" && sub_iter2 <= 
             20) {
        lambdaubold <- lambdaub
        lambdaub <- (lambdaubnew + lambdaubold)/2
        sub_iter2 <- sub_iter2 + 1
        lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                   lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                   tol = tol))
      }
      if (sub_iter2 >= 21) {
        lambdaub <- lastworking_lambdaub
        lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                   lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                   tol = tol))
        break
      }
    }
  } else if (max(lambda$lambda)/max(lambda$lambdaub) >= 1 - tol) {
    lastworking_lambdaub <- lambdaub
    while (max(lambda)/lambdaub >= 1 - tol) {
      lambdaubold <- lambdaub
      lambdaub <- lambdaubnew <- 1.2 * lambdaub
      lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                 lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                 tol = tol))
      sub_iter1 <- 1
      while (class(lambda) == "try-error" && sub_iter1 <= 
             20) {
        lambdanew <- lambdaub <- (lambdanew + lambdaold)/2
        sub_iter1 <- sub_iter1 + 1
        lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                   lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                   tol = tol))
      }
      if (sub_iter1 >= 21) {
        lambdaub <- lastworking_lambdaub
        lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                   lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                   tol = tol))
        break
      }
    }
  }
  lambdaub <- lambda$lambda
  lambda <- lambda$lambda
  
  param <- c(beta, lambda, nu, gam)
  log.Z <- logZ(log(lambda), nu, summax = summax)
  ll_new <- sum(y*log(lambda) - nu*lgamma(y+1) - log.Z)
  
  iter = iter + 1
  
  
  while ((abs((ll_new-ll_old)/ll_new)>1e-6)&& (all(abs(gam - gam_old) > 1e-06 ))&& iter <= 100){   
    ll_old <- ll_new
    nu_old <- nu
    gam_old <- gam
    log.Z <- logZ(log(lambda), nu_old, summax = summax)
    Aterm <- (comp_mean_ylogfactorialy(lambda,nu_old, log.Z, summax)-                  
                mu*comp_mean_logfactorialy(lambda,nu_old, log.Z, summax)) 
    update_score <- matrix(0, ncol = q)
    update_info_matrix <- matrix(0, ncol=t, nrow=t)
    for (k in 1:n){
      update_score <- update_score+ sum(Aterm[k] * (y[k] - mu[k])/comp_variances(lambda[k], nu_old[k], log.Z[k], summax)
                                        - (lgamma(y[k] + 1) - comp_mean_logfactorialy(lambda[k], nu_old[k], log.Z[k], summax)))*(nu_old[k]%*%S[k,]) 
      update_info_matrix <-  update_info_matrix + sum(-Aterm[k]^2/comp_variances(lambda[k], nu_old[k], log.Z[k], summax)
                                                      +  comp_variances_logfactorialy(lambda[k], nu_old[k], log.Z[k], summax))*(nu_old[k]^2*S[k,]%*%t(S[k,]))}
    gam <- gam_old + fsscale*update_score%*%solve(update_info_matrix)
    eta1 <- t(S%*%t(gam))[1,]
    nu <- exp(eta1)
    
    
    lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                               lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                               tol = tol, summax = summax))
    
    if (class(lambda) == "try-error") {
      while (class(lambda) == "try-error") {
        lambdaubold <- lambdaub
        lambdaub <- 0.8 * lambdaub
        lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                   lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                   tol = tol))
      }
      lastworking_lambdaub <- lambdaub
      sub_iter1 <- 1
      while (max(lambda)/lambdaub >= 1 - tol && sub_iter1 <= 
             20) {
        lambdaubnew <- lambdaub
        lambdaub <- (lambdaubnew + lambdaubold)/2
        sub_iter1 <- sub_iter1 + 1
        lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                   lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                   tol = tol))
        sub_iter2 <- 1
        while (class(lambda) == "try-error" && 
               sub_iter2 <= 20) {
          lambdaubold <- lambdaub
          lambdaub <- (lambdaubnew + lambdaubold)/2
          sub_iter2 <- sub_iter2 + 1
          lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                     lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                     tol = tol))
        }
        if (sub_iter2 >= 21) {
          lambdaub <- lastworking_lambdaub
          lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                     lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                     tol = tol))
          break
        }
      }
    } else if (max(lambda$lambda)/max(lambda$lambdaub) >= 1 - tol) {
      lastworking_lambdaub <- lambdaub
      while (max(lambda)/lambdaub >= 1 - tol) {
        lambdaubold <- lambdaub
        lambdaub <- lambdaubnew <- 1.2 * lambdaub
        lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                   lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                   tol = tol))
        sub_iter1 <- 1
        while (class(lambda) == "try-error" && 
               sub_iter1 <= 20) {
          lambdanew <- lambdaub <- (lambdanew + lambdaold)/2
          sub_iter1 <- sub_iter1 + 1
          lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                     lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                     tol = tol))
        }
        if (sub_iter1 >= 21) {
          lambdaub <- lastworking_lambdaub
          lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                     lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                     tol = tol))
          break
        }
      }
    }
    lambdaub <- lambda$lambdaub
    lambda <- lambda$lambda
    
    param <- c(beta, lambda, nu, gam)
    log.Z <- logZ(log(lambda), nu, summax = summax)
    ll_new <- sum(y*log(lambda) - nu*lgamma(y+1) - log.Z)
    halfstep <- 0
    while (ll_new < ll_old && halfstep <20){ 
      halfstep <- halfstep + 1
      gam <- (gam+gam_old)/2       
      eta1 <- t(S%*%t(gam))[1,]
      nu <- exp(eta1)   
      fsscale <- max(fsscale/2,1)
      lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                 lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                 tol = tol, summax = summax))
      
      if (class(lambda) == "try-error") {
        while (class(lambda) == "try-error") {
          lambdaubold <- lambdaub
          lambdaub <- 0.8 * lambdaub
          lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                     lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                     tol = tol))
        }
        lastworking_lambdaub <- lambdaub
        sub_iter1 <- 1
        while (max(lambda)/lambdaub >= 1 - tol && sub_iter1 <= 
               20) {
          lambdaubnew <- lambdaub
          lambdaub <- (lambdaubnew + lambdaubold)/2
          sub_iter1 <- sub_iter1 + 1
          lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                     lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                     tol = tol))
          sub_iter2 <- 1
          while (class(lambda) == "try-error" && 
                 sub_iter2 <= 20) {
            lambdaubold <- lambdaub
            lambdaub <- (lambdaubnew + lambdaubold)/2
            sub_iter2 <- sub_iter2 + 1
            lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                       lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                       tol = tol))
          }
          if (sub_iter2 >= 21) {
            lambdaub <- lastworking_lambdaub
            lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                       lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                       tol = tol))
            break
          }
        }
      } else if (max(lambda$lambda)/max(lambda$lambdaub) >= 1 - tol) {
        lastworking_lambdaub <- lambdaub
        while (max(lambda)/lambdaub >= 1 - tol) {
          lambdaubold <- lambdaub
          lambdaub <- lambdaubnew <- 1.2 * lambdaub
          lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                     lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                     tol = tol))
          sub_iter1 <- 1
          while (class(lambda) == "try-error" && 
                 sub_iter1 <= 20) {
            lambdanew <- lambdaub <- (lambdanew + lambdaold)/2
            sub_iter1 <- sub_iter1 + 1
            lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                       lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                       tol = tol))
          }
          if (sub_iter1 >= 21) {
            lambdaub <- lastworking_lambdaub
            lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                       lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                       tol = tol))
            break
          }
        }
      }
      lambdaub <- lambda$lambdaub
      lambda <- lambda$lambda
   
      param <- c(beta, lambda, nu, gam)
      log.Z <- logZ(log(lambda), nu, summax = summax)
      ll_new <- sum(y*log(lambda) - nu*lgamma(y+1) - log.Z)
    }                     
    
    iter = iter + 1
  }
  obj <- list()
  obj$param <- param
  obj$maxl <- ll_new
  obj$fsscale <- fsscale
  obj$iter <- iter
  obj$lambdaub <- lambdaub
  obj$update_score<- update_score                              
  obj$update_info_matrix<- update_info_matrix                             
  options(warn=0)
  return(obj)
}


