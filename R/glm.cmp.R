#' Fit a Mean Parametrized Conway-Maxwell Poisson Generalized Linear Model
#' 
#' The function \code{glm.cmp} is used to fit a mean parametrized Conway-Maxwell Poisson
#' generalized linear model with a log-link by using Fisher Scoring iteration. 
#' 
#' @usage 
#' glm.cmp(formula, data, offset = NULL,
#'    lambdalb = 1e-10, lambdaub = 1299, maxlambdaiter = 1e3, tol = 1e-6)
#' @param formula an object of class 'formula': a symblic desciption of the model to be 
#' fitted. 
#' @param data an optional data frame containing the variables in the model
#' @param offset this can be used to specify an *a priori* known component to be included 
#' in the linear predictor during fitting. This should be \code{NULL} or a numeric vector 
#' of length equal to the number of cases.  
#' @param lambdalb,lambdaub numeric: the lower and upper end points for the interval to be
#' searched for lambda(s). 
#' @param maxlambdaiter numeric: the maximum number of iterations allowed to solve 
#' for lambda(s).
#' @param tol numeric: the convergence threshold. A lambda is said to satisfy the 
#' mean constraint if the absolute difference between the calculated mean and a fitted
#' values is less than tol.
#' 
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
#' The generic assessor functions \code{coefficients} (i.e., \code{\link{coefficients}}), 
#' \code{logLik} (i.e., \code{\link{logLik.cmp}}) 
#' \code{fitted.values} (i.e., \code{\link{fitted.values}}), 
#' \code{nobs} (i.e., \code{\link{nobs.cmp}}), 
#' \code{AIC} (i.e., \code{\link{AIC.cmp}}) and 
#' \code{residuals} (i.e., \code{\link{residuals.cmp}}) 
#' can be used to extract various useful features of the value
#' returned by \code{glm.cmp}.
#' 
#' An object class 'glm.cmp' is a list containing at least the following components:
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
#' \item{terms}{the \code{\link{terms}} object used}
#' \item{data}{the \code{data} argument}
#' \item{offset}{the \code{offset} vector used}
#' 
#' @references 
#' Fung, T., Alwan, A., Wishart, J. and Huang, A. (2018). The \code{mpcmp} package for 
#' Mean-parametrized Conway-Maxwell Poisson Regression. 
#' 
#' Huang, A. (2017). Mean-parametrized Conway–Maxwell–Poisson regression models for 
#' dispersed counts. \emph{Statistical Modelling} \bold{17}, 359--380.
#'   
#' @seealso 
#' \code{\link{summary.cmp}}, \code{\link{plot.cmp}}, \code{\link{fitted.cmp}} 
#' and \code{\link{residuals.cmp}}.
#' 
#' @examples
#' ### Huang (2017) Page 368--370: Overdispersed Attendance data
#' data(attendance)
#' M.attendance = glm.cmp(daysabs~ gender+math+prog, data=attendance)
#' M.attehdance
#' summary(M.attendance)
#' 
#' ### Huang (2017) Page 371--372: Underdispersed Takeover Bids data
#' data(takeoverbids)
#' M.bids = glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#' + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#' M.bids
#' summary(M.bids)


glm.cmp <- function(formula, data, offset = NULL,
                    lambdalb = 1e-10, lambdaub = 1299, maxlambdaiter = 1e3, tol = 1e-6){
  call <- match.call()
  if (is.null(formula)) {
    stop("formula must be specified (can not be NULL)")
  }
  if (lambdalb>=lambdaub) {
    stop("lower bound for the search of lambda must be smaller than the upper bound")
  }
  if (missing(data)){
    data <- environment(formula)
  }
  mf <- stats::model.frame(formula, data=data)
  mt <- attr(mf, "terms")
  y <- model.response(mf)
  X <- model.matrix(formula,mf)
  if (is.null(offset)){
    offset.cmp = rep(0,length(y))
  } else {
    offset.cmp = offset
  }
  # use poisson glm to generate initial values for betas
  M0 <- glm(y~-1+X+offset(offset.cmp), family=poisson())
  offset <- M0$offset
  n <- length(y) # sample size
  q <- ncol(X)  # number of covariates
  #starting values for optimization
  beta0 <- coef(M0)
  nu0 <- 1
  nu_lb <- 1e-10
  lambda0 <- fitted.values(M0)
  param <- c(beta0,lambda0,nu0)
  ll_old <- comp_mu_loglik(param = param, y=y, xx= X, offset=offset)
  param_obj<- getnu(param = param, y=y, xx= X, offset = offset, llstart = ll_old, fsscale=32)
  ll_new <- param_obj$maxl
  param <- param_obj$param
  fsscale <- param_obj$fsscale
  iter <- 0
  while (abs((ll_new-ll_old)/ll_new)>1e-6 && iter<=100){
    iter <- iter +1
    ll_old <- ll_new
    paramold <- param
    betaold <- param[1:q]
    etaold <- t(X%*%betaold)[1,]
    muold <-  exp(etaold+offset)
    lambdaold <- param[(q+1):(q+n)]
    nuold <- param[q+n+1]
    W <- diag(muold^2/comp_variances(lambdaold,nuold))
    z <- etaold + (y-muold)/muold
    beta <- solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%z
    Aterm <- (comp_mean_ylogfactorialy(lambdaold,nuold)-
                muold*comp_mean_logfactorialy(lambdaold,nuold))
    update_score <-
      sum(Aterm*(y-muold)/comp_variances(lambdaold,nuold)
          -(lgamma(y+1)-comp_mean_logfactorialy(lambdaold,nuold)))
    update_info_matrix <-
      sum(Aterm/comp_variances(lambdaold,nuold)+
            comp_variances_logfactorialy(lambdaold,nuold))
    nu <- nuold + update_score/update_info_matrix
    while (nu < nu_lb){
      nu <- (nu+nuold)/2
    }
    eta <- t(X%*%beta)[1,]
    mu <- exp(eta+offset)
    lambda <- comp_lambdas(mu,nu)
    param <- c(beta, lambda, nu)
    ll_new <- comp_mu_loglik(param = param, y=y, xx= X, offset= offset)
    halfstep <- 0
    while (ll_new < ll_old && halfstep <= 20){
      halfstep <- halfstep + 1
      beta <- (beta+betaold)/2
      nu <- (nu+nuold)/2
      eta <- t(X%*%beta)[1,]
      mu <- exp(eta+offset)
      lambda <- comp_lambdas(mu,nu)
      param <- c(beta, lambda, nu)
      ll_new <- comp_mu_loglik(param = param, y=y, xx= X, offset= offset)
    }
  }
  maxl = ll_new # maximum loglikelihood achieved
  beta = param[1:q]  # estimated regression coefficients beta
  lambda = param[(q+1):(q+n)]  # estimated rates (not generally useful)
  nu = param[q+n+1] # estimate dispersion
  precision_beta = 0
  eta <- t(X%*%beta)[1,]
  if (is.null(offset)){
    fitted <- exp(eta)
  } else {
    fitted <- exp(eta+offset)
  }
  for (i in 1:n){
    precision_beta = precision_beta + fitted[i]^2*X[i,]%*%t(X[i,])/comp_variances(lambda[i], nu)
  }
  variance_beta <- solve(precision_beta)
  se_beta <- as.vector(sqrt(diag(variance_beta)))
  Xtilde <- diag(fitted/sqrt(comp_variances(lambda,nu)))%*%as.matrix(X)
  h <- diag(Xtilde%*%solve(t(Xtilde)%*%Xtilde)%*%t(Xtilde))
  df.residuals <- length(y)-length(beta)
  if (df.residuals > 0){
    indsat.deviance = dcomp(y, mu = y, nu = nu, log.p=TRUE, lambdalb = lambdalb,
                            lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, tol = tol)
    indred.deviance = 2*(indsat.deviance - dcomp(y, nu=nu, lambda= lambda,
                                                 log.p=TRUE))
    d.res = sign(y-fitted)*sqrt(abs(indred.deviance))
  } else { d.res = rep(0,length(y))
  }
  out <- list()
  out$call <- call
  out$formula <- formula
  out$y <- y
  out$x <- X
  out$data <- data
  out$nobs <- n
  out$iter <- iter
  out$coefficients <- beta
  out$rank <- length(beta)
  out$lambda <- lambda
  out$offset <- offset
  out$nu <- nu
  out$terms <- mt
  out$model <- mf
  out$linear.predictors <- eta
  out$maxl <- maxl
  out$fitted.values <- fitted
  out$residuals <- y - fitted
  out$leverage <- h
  out$d.res <- d.res
  out$variance_beta <- variance_beta
  out$se_beta <- se_beta
  out$df.residuals <- df.residuals
  out$df.null <- n-1
  out$null.deviance <- 2*(sum(indsat.deviance) -
                            sum(dcomp(y, lambda =
                                        comp_lambdas(mean(y), nu,
                                                     lambdalb = lambdalb,
                                                     lambdaub = lambdaub,
                                                     maxlambdaiter = maxlambdaiter,
                                                     tol = tol),
                                      nu = nu,log.p = TRUE)))
  out$residuals.deviance <- 2*(sum(indsat.deviance) -
                                 sum(dcomp(y, lambda = lambda,
                                           nu = nu,log.p = TRUE)))
  names(out$coefficients) = labels(X)[[2]]
  class(out) <- "cmp"
  return(out)
}
