#' Parameter Generator for nu
#' 
#' A function that use the arguments of a \code{glm.cmp} call to generate a better initial 
#' \code{nu} estimate. 
#' 
#' @param param numeric vector: the model coefficients & the current value of \code{nu}. 
#' It is assumed that \code{nu} is in the last position of \code{param}.
#' @param y numeric vector: response variable
#' @param xx numeric matrix: the explanatory variables
#' @param offset numeric vector: a vector of lenght equal to the number of cases
#' @param llstart numeric: current log-likelihood value
#' @param fsscale numeric: a scaling factor (generally >1) for the 
#' relaxed fisher scoring algorithm
#' @return 
#' List containing the following:
#' \item{param}{the model coefficients & the updated \code{nu}}
#' \item{maxl}{the updated log-likelihood}
#' \item{fsscale}{the final scaling factor used}
#' 
getnu <- function(param, y, xx , offset, llstart, fsscale = 1){
  n <- length(y) # sample size
  q <- ncol(xx)  # number of covariates
  nu_lb <- 1e-10
  iter <- 1
  beta <- param[1:q]
  mu <-  exp(t(xx%*%beta)[1,])
  lambda <- param[(q+1):(q+n)]
  nu_old <- param[q+n+1]
  ll_old <- llstart
  Aterm <- (comp_mean_ylogfactorialy(lambda,nu_old)-
              mu*comp_mean_logfactorialy(lambda,nu_old))
  update_score <-
    sum(Aterm*(y-mu)/comp_variances(lambda,nu_old)
        -(lgamma(y+1)-comp_mean_logfactorialy(lambda,nu_old)))
  update_info_matrix <-
    sum(Aterm/comp_variances(lambda,nu_old)+
          comp_variances_logfactorialy(lambda,nu_old))
  nu <- nu_old + fsscale*update_score/update_info_matrix
  while (nu < nu_lb){
    nu <- (nu+nu_old)/2
  }
  lambda <- comp_lambdas(mu,nu)
  param <- c(beta, lambda, nu)
  ll_new <- comp_mu_loglik(param = param, y=y, xx= xx, offset= offset)
  while (ll_new < ll_old){
    nu <- (nu+nu_old)/2
    fsscale <- max(fsscale/2,1)
    lambda <- comp_lambdas(mu,nu)
    param <- c(beta, lambda, nu)
    ll_new <- comp_mu_loglik(param = param, y=y, xx= xx, offset= offset)
  }
  iter = iter+1
  while ((abs((ll_new-ll_old)/ll_new)>1e-6) && abs(nu-nu_old)>1e-6 && iter <=100){
    ll_old <- ll_new
    nu_old <- nu
    Aterm <- (comp_mean_ylogfactorialy(lambda,nu_old)-
                mu*comp_mean_logfactorialy(lambda,nu_old))
    update_score <-
      sum(Aterm*(y-mu)/comp_variances(lambda,nu_old)
          -(lgamma(y+1)-comp_mean_logfactorialy(lambda,nu_old)))
    update_info_matrix <-
      sum(Aterm/comp_variances(lambda,nu_old)+
            comp_variances_logfactorialy(lambda,nu_old))
    nu <- nu_old + fsscale*update_score/update_info_matrix
    while (nu < nu_lb){
      nu <- (nu + nu_old)/2
    }
    if (abs(nu-nu_old)<1e-6){
      break
    }
    lambda <- comp_lambdas(mu,nu)
    param <- c(beta, lambda, nu)
    ll_new <- comp_mu_loglik(param = param, y=y, xx= xx, offset= offset)
    while (ll_new < ll_old){
      nu <- (nu+nu_old)/2
      fsscale <- max(fsscale/2,1)
      lambda <- comp_lambdas(mu,nu)
      param <- c(beta, lambda, nu)
      ll_new <- comp_mu_loglik(param = param, y=y, xx= xx, offset= offset)
    }
    iter = iter+1
  }
  obj <- list()
  obj$param = param
  obj$maxl = ll_new
  obj$fsscale = fsscale
  return(obj)
}
