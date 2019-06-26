#' Fit a Mean Parametrized Conway-Maxwell Poisson Generalized Linear Model
#' 
#' The function \code{glm.cmp} is used to fit a mean parametrized Conway-Maxwell Poisson
#' generalized linear model with a log-link by using Fisher Scoring iteration. 
#' 
#' @usage 
#' glm.cmp(formula, data, offset = NULL, subset, na.action, betastart = NULL,
#'    lambdalb = 1e-10, lambdaub = 1000, maxlambdaiter = 1e3, tol = 1e-6, 
#'    contrasts = NULL)
#' @param formula an object of class 'formula': a symblic desciption of the model to be 
#' fitted. 
#' @param data an optional data frame containing the variables in the model
#' @param offset this can be used to specify an *a priori* known component to be included 
#' in the linear predictor during fitting. This should be \code{NULL} or a numeric vector 
#' of length equal to the number of cases.  
#' @param subset an optional vector specifying a subset of observations to be used in the 
#' fitting process.
#' @param na.action a function which indicates what should happen when the data contain 
#' NAs. The default is set by the na.action setting of options, and is na.fail if that 
#' is unset. The ‘factory-fresh’ default is na.omit. Another possible value is NULL, 
#' no action. Value na.exclude can be useful.
#' @param betastart starting values for the parameters in the linear predictor for mu.
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
#' Y_i ~ CMP(mu_i, nu), 
#'           
#' where  
#'    
#' E(Y_i) = mu_i = exp(x_i^T beta),
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

glm.cmp <- function(formula, data, offset = NULL,
                    subset, na.action, betastart = NULL, 
                    lambdalb = 1e-10, lambdaub = 1000, maxlambdaiter = 1e3, tol = 1e-6,
                    contrasts = NULL){
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", 
               "offset"), names(mf), 0L)
  if (is.null(formula)) {
    stop("formula must be specified (can not be NULL)")
  }
  if (lambdalb>=lambdaub) {
    stop("lower bound for the search of lambda must be smaller than the upper bound")
  }
  if (missing(data)){
    data <- environment(formula)
  }
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf)
  X <- model.matrix(formula,mf,contrasts)
  if (is.null(offset)){
    offset.cmp = rep(0,length(y))
  } else {
    offset.cmp = offset
  }
  # use poisson glm to generate initial values for betas
  M0 <- stats::glm(y~-1+X+offset(offset.cmp), start = betastart, family=stats::poisson())
  offset <- M0$offset
  n <- length(y) # sample size
  q <- ncol(X)  # number of covariates
  #starting values for optimization
  beta0 <- stats::coef(M0)
  mu0 <- M0$fitted.values
  nu_lb <- 1e-10
  summax <- ceiling(max(c(max(y)+20*sqrt(var(y)),100)))
  #summax <- 100
  nu0 <- exp(optimize(f= comp_mu_loglik_log_nu_only, 
                      interval = c(max(-log(1+2*mu0)), 5), mu=mu0, y=y, 
                      summax = summax)$minimum)
  lambda0 <- (mu0+(nu0-1)/(2*nu0))^(nu0)
  summax <- ceiling(max(c(mu0+20*sqrt(mu0/nu0),100)))
  #lambdaub <- min(lambdaub, 2*max(lambda0))
  lambda.ok <- comp_lambdas(mu0, nu0, lambdalb = lambdalb, 
                            #lambdaub = lambdaub,
                            lambdaub = min(lambdaub,2*max(lambda0)), 
                            maxlambdaiter = maxlambdaiter, tol = tol, summax = summax, 
                            lambdaint = lambda0)
  lambda0 <- lambda.ok$lambda
  lambdaub <- lambda.ok$lambdaub
  param <- c(beta0,lambda0,nu0)
  ll_old <- comp_mu_loglik(param = param, y=y, xx= X, offset=offset, summax = summax)
  param_obj<- getnu(param = param, y=y, xx= X, offset = offset, llstart = ll_old, 
                    fsscale = 1, lambdalb = lambdalb, lambdaub = lambdaub, 
                    maxlambdaiter = maxlambdaiter, tol = tol,summax = summax)
  lambdaub <- param_obj$lambdaub 
  ll_new <- param_obj$maxl
  param <- param_obj$param
  fsscale <- param_obj$fsscale
  iter <- param_obj$iter
  while (abs((ll_new-ll_old)/ll_new) > tol && iter <= 100){
    iter <- iter +1
    ll_old <- ll_new
    paramold <- param
    betaold <- param[1:q]
    etaold <- t(X%*%betaold)[1,]
    muold <-  exp(etaold+offset)
    lambdaold <- param[(q+1):(q+n)]
    nuold <- param[q+n+1]
    log.Z <- logZ(log(lambdaold), nuold, summax = summax)
    ylogfactorialy <- comp_mean_ylogfactorialy(lambdaold, nuold, log.Z, summax)
    logfactorialy <- comp_mean_logfactorialy(lambdaold, nuold, log.Z, summax)
    variances <- comp_variances(lambdaold, nuold, log.Z, summax)
    variances_logfactorialy <- 
      comp_variances_logfactorialy(lambdaold, nuold, log.Z, summax)
    W <- diag(muold^2/variances)
    z <- etaold + (y-muold)/muold
    beta <- solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%z
    eta <- t(X%*%beta)[1,]
    mu <- exp(eta+offset)
    Aterm <- (ylogfactorialy- mu*logfactorialy)
    update_score <- sum(Aterm*(y-mu)/variances -(lgamma(y+1)-logfactorialy))
    update_info_matrix <- sum(-Aterm^2/variances+ variances_logfactorialy)
    if (update_info_matrix < 0){
      update_info_matrix <- sum(variances_logfactorialy)
    }
    nu <- nuold + update_score/update_info_matrix
    while (nu < nu_lb){
      nu <- (nu+nuold)/2
    }
    lambda.ok <- comp_lambdas(mu,nu, lambdalb = lambdalb,
                              #lambdaub = lambdaub,
                              lambdaub = min(lambdaub,2*max(lambdaold)), 
                              maxlambdaiter = maxlambdaiter, tol = tol,
                              lambdaint = lambdaold, summax = summax)
    lambda <- lambda.ok$lambda
    lambdaub <- lambda.ok$lambdaub
    param <- c(beta, lambda, nu)
    ll_new <- comp_mu_loglik(param = param, y=y, xx= X, offset= offset, summax = summax)
    halfstep <- 0
    while (ll_new < ll_old && halfstep <= 20 && abs((ll_new-ll_old)/ll_new)>tol){
      halfstep <- halfstep + 1
      beta <- (beta+betaold)/2
      nu <- (nu+nuold)/2
      eta <- t(X%*%beta)[1,]
      mu <- exp(eta+offset)
      lambdaold <- lambda
      lambda.ok <- comp_lambdas(mu,nu, lambdalb = lambdalb,
                                #lambdaub = lambdaub,
                                lambdaub = min(lambdaub,2*max(lambdaold)), 
                                maxlambdaiter = maxlambdaiter, tol = tol,
                                lambdaint = lambda, summax = summax)
      lambda <- lambda.ok$lambda
      lambdaub <- lambda.ok$lambdaub
      param <- c(beta, lambda, nu)
      ll_new <- comp_mu_loglik(param = param, y=y, xx= X, offset= offset, summax)
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
  log.Z <- logZ(log(lambda), nu, summax = summax)
  variances <- comp_variances(lambda, nu, log.Z = log.Z, summax = summax)
  for (i in 1:n){
    precision_beta = precision_beta + fitted[i]^2*X[i,]%*%t(X[i,])/variances[i]
  }
  variance_beta <- solve(precision_beta)
  se_beta <- as.vector(sqrt(diag(variance_beta)))
  Xtilde <- diag(fitted/sqrt(variances))%*%as.matrix(X)
  h <- diag(Xtilde%*%solve(t(Xtilde)%*%Xtilde)%*%t(Xtilde))
  df.residuals <- length(y)-length(beta)
  if (df.residuals > 0){
    indsat.deviance = dcomp(y, mu = y, nu = nu, log.p=TRUE, lambdalb = min(lambdalb),
                            lambdaub = lambdaub, maxlambdaiter = maxlambdaiter, 
                            tol = tol, summax =summax)
    indred.deviance = 2*(indsat.deviance - dcomp(y, nu=nu, lambda= lambda,
                                                 log.p=TRUE, summax=summax))
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
  out$log.Z <- log.Z
  out$summax <- summax 
  out$offset <- offset
  out$nu <- nu
  out$terms <- mt
  out$model <- mf
  out$contrasts <- attr(X, "contrasts")
  out$na.action <- attr(mf, "na.action")
  out$lambdaub <- lambdaub
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