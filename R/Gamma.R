
#' Fit a Mean Parametrized Conway-Maxwell Poisson Generalized Linear Model with varying dispersion
#' 
#' The function \code{glm.cmp.gamma} is used to fit a mean parametrized Conway-Maxwell Poisson
#' generalized linear model with a log-link by using Fisher Scoring iteration. 
#' 
#' @usage 
#' glm.cmp.gamma(formula, dformula, data, offset = NULL,  betastart = NULL,
#'    ambdalb = 1e-10, lambdaub = 1900, maxlambdaiter = 1000, 
#'    contrasts = NULL)
#' @param formula an object of class 'formula': a symblic desciption of the model to be 
#' fitted. 
#' @param dformula an object of class 'dformula': a symblic desciption of the model for dispersion.
#' @param data an optional data frame containing the variables in the model.
#' @param offset this can be used to specify an *a priori* known component to be included 
#' in the linear predictor during fitting. This should be \code{NULL} or a numeric vector 
#' of length equal to the number of cases.  
#' @param betastart starting values for the parameters in the linear predictor for mu.
#' @param gammastart starting values for the parameters in the linear predictor for nu.
#' @param lambdalb,lambdaub numeric: the lower and upper end points for the interval to be
#' searched for lambda(s). The default value for lambdaub should be sufficient for small to
#' moderate size nu. If nu is large and required a larger \code{lambdaub}, the algorithm
#' will scale up \code{lambdaub} accordingly.  
#' @param maxlambdaiter numeric: the maximum number of iterations allowed to solve 
#' for lambda(s).
#' @param tol numeric: the convergence threshold. A lambda is said to satisfy the 
#' mean constraint if the absolute difference between the calculated mean and a fitted
#' values is less than tol.
#' @param contrasts an optional list. See the contrasts.arg of model.matrix.default.
#' @export
#' @import stats
#' @details 
#' Fit a mean-parametrizied COM-Poisson regression using maximum likelihood estimation 
#' via an iterative Fisher Scoring algorithm. 
#' 
#' The COM-Poisson regression model is
#' 
#' Y_i ~ CMP(mu_i, nu_i), 
#'           
#' where  
#'    
#' E(Y_i) = mu_i = exp(x_i^T beta),
#'  nu_i = exp(x_i^T gamma),
#'       
#' and \emph{nu > 0} is the dispersion parameter. 
#' 
#' The fitted COM-Poisson distribution is over- or under-dispersed 
#' if \emph{nu < 1} and \emph{nu > 1} respectively.
#' @return 
#' A fitted model object of class \code{cmp} similar to one obtained from \code{glm} 
#' or \code{glm.nb}.
#' 
#' The function \code{summary} (i.e., \code{\link{summary.cmp}}) can be used to obtain 
#' and print a summary of the results. 
#' 
#' The function \code{plot} (i.e., \code{\link{plot.cmp}}) can be used to produce a range 
#' of diagnostic plots. 
#' 
#' The generic assessor functions \code{coef} (i.e., \code{\link{coef.cmp}}), 
#' \code{logLik} (i.e., \code{\link{logLik.cmp}}) 
#' \code{fitted} (i.e., \code{\link{fitted.cmp}}), 
#' \code{nobs} (i.e., \code{\link{nobs.cmp}}), 
#' \code{AIC} (i.e., \code{\link{AIC.cmp}}) and 
#' \code{residuals} (i.e., \code{\link{residuals.cmp}}) 
#' can be used to extract various useful features of the value
#' returned by \code{glm.cmp.gamma}.
#' 
#' An object class 'glm.cmp.gamma' is a list containing at least the following components:
#'
#' \item{coefficients}{a named vector of coefficients}
#' \item{se_beta}{approximate standard errors (using observed rather than expected 
#' information) for coefficients}
#' \item{residuals}{the \emph{response} residuals (i.e., observed-fitted)}
#' \item{fitted.values}{the fitted mean values}
#' \item{rank}{the numeric rank of the fitted linear model}
#' \item{linear.predictors}{the linear fit on log scale}
#' \item{df.residuals}{the residuals degrees of freedom}
#' \item{df.null}{the residual degrees of freedom for the null model}
#' \item{null.deviance}{The deviance for the null model. 
#' The null model will include only the intercept.}
#' \item{y}{the \code{y} vector used.}
#' \item{x}{the model matrix}
#' \item{model}{the model frame}
#' \item{call}{the matched call}
#' \item{formula}{the formula supplied}
#' \item{terms}{the \code{terms} object used}
#' \item{data}{the \code{data} argument}
#' \item{offset}{the \code{offset} vector used}
#' \item{lambdaub}{the final \code{lambdaub} used}
#' 
#' @references 
#' Fung, T., Alwan, A., Wishart, J. and Huang, A. (2019). \code{mpcmp}: Mean-parametrized
#' Conway-Maxwell Poisson Regression. R package version 0.2.0.
#' 
#' Huang, A. (2017). Mean-parametrized Conway-Maxwell-Poisson regression models for 
#' dispersed counts. \emph{Statistical Modelling} \bold{17}, 359--380.
#'   
#' @seealso 
#' \code{\link{summary.cmp}}, \code{\link{plot.cmp}}, \code{\link{fitted.cmp}} 
#' and \code{\link{residuals.cmp}}.
#' @examples 
#' ### Huang (2017) Page 368--370: Overdispersed Attendance data
#' data(attendance)
#' M.attendance <- glm.cmp(daysabs~ gender+math+prog, data=attendance)
#' M.attendance
#' summary(M.attendance)
#' plot(M.attendance)
#' 
#' ### Barbour & Brown (1974): Overdispersed Fish data
#' data(fish)
#' M.fish <- glm.cmp(species~ 1+log(area), data=fish)
#' M.fish
#' summary(M.fish)
#' 
#' ### Huang (2017) Page 371--372: Underdispersed Takeover Bids data
#' data(takeoverbids)
#' M.bids <- glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#'     + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#' M.bids
#' summary(M.bids)
#' par(mfrow=c(2,2))
#' plot(M.bids)
#' 
#' ### Huang (2017) Page 373--375: Underdispersed Cotton bolls data
#' ### Model fitting for predictor V 
#' \donttest{
#' data(cottonbolls)
#' M.bolls <- glm.cmp(nc~ 1+stages:def+stages:def2, data= cottonbolls)
#' M.bolls
#' summary(M.bolls)
#' }
glm.cmp.gamma<- function (formula, dformula, data, offset = NULL, betastart = NULL,   #dataS,
                         lambdalb = 1e-10, lambdaub = 1900, maxlambdaiter = 1000, 
                         tol = 1e-06, contrasts = NULL) 
{
  
  if (missing(data)){ 
    data <- environment(formula)
  }
  frame <- model.frame(formula, data)
  terms <- attr(frame, "terms")
  X <- model.matrix(terms, frame)
  S <- model.matrix(terms, frame)
  y <- model.response(frame)
  X <- model.matrix(formula,data)                   
  S <- model.matrix(dformula, data)                 
  q <- ncol(X)
  t <- ncol(S)
  if (dformula == ~1) 
    colnames(S) <- "nu"
  if (is.null(offset)) {
    offset.cmp = rep(0, length(y))
  }else {
    offset.cmp = offset
  }
  M0 <- stats::glm(y ~ -1 + X + offset(offset.cmp), start = betastart, family = stats::poisson())
  offset <- M0$offset
  nu_lb <- 1e-10
  summax <- 300
                   
  n <- length(y)
  q <- ncol(X)
  t <- ncol(S)
  beta0 <- stats::coef(M0)                           
  gam0 <-  rep(0,ncol(S))
  eta10 <- t(S %*% gam0)[1, ]                                                          
  nu0 <- exp(eta10 + offset)  
  lambda0 <- stats::fitted(M0)
  param <- c(beta0, lambda0, nu0, gam0)                                             
  log.Z <- logZ(log(lambda0), nu0, summax = summax)
  ll_old <- sum(y*log(lambda0) - nu0*lgamma(y+1) - log.Z)
  param_obj<- getnu.gamma(param, y, X , S = S, offset, llstart = ll_old, fsscale = 1, 
                          lambdalb = 1e-10, lambdaub = 1900, maxlambdaiter = 1e3, tol = 1e-6)  
  lambdaub <- param_obj$lambdaub
  ll_new <- param_obj$maxl
  param <- param_obj$param
  fsscale <- 1 
  iter <- 0
 
  while (abs((ll_new - ll_old)/ll_new) > tol && iter <= 100) {
    iter <- iter + 1
    ll_old <- ll_new
    paramold <- param
    betaold <- param[1:q]
    etaold <- t(X %*% betaold)[1, ]
    muold <- exp(etaold + offset)
    lambdaold <- param[(q + 1):(q + n)]
    nuold <- param[(q+n+1):(q+2*n)]                                      
    gamold <- param[(q+2*n+1):(q+2*n+t)]
    log.Z <- logZ(log(lambdaold), nuold, summax = summax)
    W <- diag(muold^2/comp_variances(lambdaold, nuold, log.Z))
    z <- etaold + (y - muold)/muold
    beta <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% z
    update_score <- matrix(0, ncol = q)
    update_info_matrix <- matrix(0, ncol=t, nrow=t)
    Aterm <- (comp_mean_ylogfactorialy(lambdaold, nuold, log.Z, summax) - 
                muold * comp_mean_logfactorialy(lambdaold, nuold, log.Z, summax))  #n*1
    for (k in 1:n){
      update_score <- update_score+ sum(Aterm[k] * (y[k] - muold[k])/comp_variances(lambdaold[k], nuold[k], log.Z[k], summax)
                   - (lgamma(y[k] + 1) - comp_mean_logfactorialy(lambdaold[k], nuold[k], log.Z[k], summax)))*(nuold[k]%*%S[k,]) 
      update_info_matrix <- update_info_matrix+  sum(-Aterm[k]^2 /comp_variances(lambdaold[k], nuold[k], log.Z[k], summax)
                     + comp_variances_logfactorialy(lambdaold[k], nuold[k], log.Z[k], summax))*(nuold[k]^2*S[k,]%*%t(S[k,]))}
    gam <- gamold + update_score%*%solve(update_info_matrix)           
    eta1 <- t(S %*% t(gam))[1, ]                                     
    nu <- exp(eta1)
    
    eta <- t(X %*% beta)[1, ]
    mu <- exp(eta + offset)
    lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                               lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                               tol = tol, summax = summax))
    
    
    
    if (class(lambda) == "try-error") {
      while (class(lambda) == "try-error") {
        lambdaubold <- lambdaub
        lambdaub <- 0.8 * lambdaub
        lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                   lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                   tol = tol, summax = summax))
        lambda <- lambda.ok$lambda
        lambdaub <-lambda.ok$lambdaub
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
                                   tol = tol, summax= summax))
        sub_iter2 <- 1
        while (class(lambda) == "try-error" && 
               sub_iter2 <= 20) {
          lambdaubold <- lambdaub
          lambdaub <- (lambdaubnew + lambdaubold)/2
          sub_iter2 <- sub_iter2 + 1
          lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                     lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                     tol = tol, summax = summax))
        }
        if (sub_iter2 >= 21) {
          lambdaub <- lastworking_lambdaub
          lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                     lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                     tol = tol, summax = summax))
          break
        }
      }
    } else if (max(lambda$lambda)/max(lambda$lambdaub) >= 1 - tol) {
      lastworking_lambdaub <- lambdaub
      while (max(lambda)/lambdaub >= 1 - tol) {
        lambdaubold <- lambdaub
        lambdaub <- lambdaubnew <- 1.2 * lambdaub
        lambda.ok <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                   lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                   tol = tol, summax = summax))
        sub_iter1 <- 1
        while (class(lambda) == "try-error" && 
               sub_iter1 <= 20) {
          lambdanew <- lambdaub <- (lambdanew + lambdaold)/2
          sub_iter1 <- sub_iter1 + 1
          lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                     lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                     tol = tol, summax = summax))
        }
        if (sub_iter1 >= 21) {
          lambdaub <- lastworking_lambdaub
          lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                     lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                     tol = tol, summax=summax))
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
    while (ll_new < ll_old && halfstep <= 20) {
      halfstep <- halfstep + 1
      beta <- (beta + betaold)/2
      nu <- (nu + nuold)/2
      eta <- t(X %*% beta)[1, ]
      mu <- exp(eta + offset)
      lambda <- try(comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                    lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                    tol = tol, summax = summax))
      if (class(lambda) == "try-error") {
        while (class(lambda) == "try-error") {
          lambdaubold <- lambdaub
          lambdaub <- 0.8 * lambdaub
          lambda <- try(mpcmp:::comp_lambdas(mu, nu, lambdalb = lambdalb, 
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
          lambda <- try(mpcmp:::comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                             lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                             tol = tol))
          sub_iter2 <- 1
          while (class(lambda) == "try-error" && 
                 sub_iter2 <= 20) {
            lambdaubold <- lambdaub
            lambdaub <- (lambdaubnew + lambdaubold)/2
            sub_iter2 <- sub_iter2 + 1
            lambda <- try(mpcmp:::comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                               lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                               tol = tol))
          }
          if (sub_iter2 >= 21) {
            lambdaub <- lastworking_lambdaub
            lambda <- try(mpcmp:::comp_lambdas(mu, nu, lambdalb = lambdalb, 
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
          lambda <- try(mpcmp:::comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                             lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                             tol = tol))
          sub_iter1 <- 1
          while (class(lambda) == "try-error" && 
                 sub_iter1 <= 20) {
            lambdanew <- lambdaub <- (lambdanew + lambdaold)/2
            sub_iter1 <- sub_iter1 + 1
            lambda <- try(mpcmp:::comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                               lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                               tol = tol))
          }
          if (sub_iter1 >= 21) {
            lambdaub <- lastworking_lambdaub
            lambda <- try(mpcmp:::comp_lambdas(mu, nu, lambdalb = lambdalb, 
                                               lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                                               tol = tol))
            break
          }
        }
      }
      lambdaub <- lambda$lambdaub
      lambda <- lambda$lambda
      param <- c(beta, lambda, nu,gam)
      log.Z <- logZ(log(lambda), nu, summax = summax)
      ll_new <- sum(y*log(lambda) - nu*lgamma(y+1) - log.Z)
    }
  }
  maxl = ll_new
  beta = param[1:q]
  lambda = param[(q + 1):(q + n)]
  gam <- param[(q+2*n+1):(q+2*n+t)]             
  nu <- param[(q+n+1):(q+2*n)]                        
  #nu = param[q + n + 1]
  precision_beta = 0
  eta <- t(X %*% beta)[1, ]
  eta1 <- t(S%*%gam)[1,]
  nu <- exp(eta1)
  if (is.null(offset)) {
    fitted <- exp(eta)
  } else {
    fitted <- exp(eta + offset)
  }
  for (i in 1:n) {
    precision_beta = precision_beta + fitted[i]^2 * X[i, ] %*% t(X[i, ])/mpcmp:::comp_variances(lambda[i], nu[i]) 
  }
  variance_beta <- solve(precision_beta)
  se_beta <- as.vector(sqrt(diag(variance_beta)))
  Xtilde <- diag(fitted/sqrt(mpcmp:::comp_variances(lambda, nu))) %*% as.matrix(X)
  h <- diag(Xtilde %*% solve(t(Xtilde) %*% Xtilde) %*% t(Xtilde))
  df.residuals <- length(y) - length(beta)
  if (df.residuals > 0) {
    indsat.deviance = mpcmp:::dcomp(y, mu = y, nu = nu, log.p = TRUE, 
                            lambdalb = lambdalb, lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                            tol = tol)
    indred.deviance = 2 * (indsat.deviance - mpcmp:::dcomp(y, nu = nu, 
                                                   lambda = lambda, log.p = TRUE))
    d.res = sign(y - fitted) * sqrt(abs(indred.deviance))
  }else {
    d.res = rep(0, length(y))
  }

  out <- list()
  out$call <- call
  out$formula <- formula
  out$dformula <- dformula
  out$y <- y
  out$x <- X
  out$S <- S
  out$data <- data                           
  out$nobs <- n
  out$iter <- iter
  out$coefficients <- beta
  out$rank <- length(beta)
  out$lambda <- lambda
  out$offset <- offset
  out$gam <- gam
  out$nu <- nu
  out$beta < beta
  out$contrasts <- attr(X, "contrasts")
  out$lambdaub <- lambdaub
  out$linear.predictors <- eta
  out$linear.predictors <- eta1
  out$maxl <- maxl
  out$fitted.values <- fitted
  out$residuals <- y - fitted
  out$leverage <- h
  out$d.res <- d.res
  out$variance_beta <- variance_beta
  out$se_beta <- se_beta
  out$df.residuals <- df.residuals
  out$df.null <- n - 1
  
  out$null.deviance <- 2*(sum(indsat.deviance) -
                            sum(dcomp(y, mu = mean(y), nu = nu, log.p = TRUE, 
                                      summax=summax, lambdalb = min(lambdalb),
                                      lambdaub = max(lambdaub))))
  out$residuals.deviance <- 2*(sum(indsat.deviance) -
                                 sum(dcomp(y, lambda = lambda, nu = nu, 
                                           log.p = TRUE, summax=summax)))
  names(out$coefficients) = labels(X)[[2]]
  class(out) <- "cmp"
  return(out)
}
