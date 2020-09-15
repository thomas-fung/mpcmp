#' Fit a Mean Parametrized Conway-Maxwell Poisson Generalized Linear 
#' Model with varying dispersion. 
#'
#' @param y the response y vector.
#' @param X the design matrix for regressing the mean
#' @param S the design matrix for regressing the dispersion
#' @param offset this can be used to specify an *a priori* known component to be included 
#' in the linear predictor for mean during fitting. This should be \code{NULL} or a numeric vector 
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
#'
#' @return
#' A fitted model object of class \code{cmp} similar to one obtained from \code{glm} or \code{glm.nb}.
#'
#' @examples
#' ## For examples see example(glm.cmp)
fit_glm_cmp_vary_nu <- function(y=y, X = X, S = S, offset = offset,
                                betastart = betastart, 
                                gammastart = gammastart,
                                lambdalb = lambdalb, lambdaub = lambdaub, 
                                maxlambdaiter = maxlambdaiter, tol = tol) {
  M0 <- stats::glm(y~-1+X+offset(offset), start = betastart, 
                   family=stats::poisson())
  offset <- M0$offset
  beta0 <- stats::coef(M0)
  n <- length(y) # sample size
  q1 <- ncol(X)  # number of covariates for mu
  q2 <- ncol(S)  # number of covariates for nu
  summax <- ceiling(max(c(max(y)+20*sqrt(var(y)),100)))
  if (!is.null(gammastart) && q2 != length(gammastart)){
    stop(paste("length of 'gammastart' should equal to", q2, 
               "\nand corresponding to initial coefs for ", 
               paste(colnames(S), collapse =", ")))
  } else if (is.null(gammastart)){
    gamma0 <- rep(0, q2)
    nu0 <- rep(1, n)
    lambda0 <- mu0 <- M0$fitted.values
  } else {
    gamma0 <- gammastart 
    nu0 <- exp(t(S%*%gamma0)[1,])
    mu0 <- M0$fitted.values
    lambda.ok <- comp_lambdas(mu0, nu0, lambdalb = lambdalb, 
                              lambdaub = min(lambdaub, max(lambdaold*2)), 
                              maxlambdaiter = maxlambdaiter, tol = tol,
                              lambdaint = lambdaold, summax = summax)
    lambda0 <- lambda.ok$lambda
    lambdaub <-lambda.ok$lambdaub
  }
  #starting values for optimization
  param <- c(beta0, lambda0, nu0, gamma0) 
  # computing log-likelihood 
  ll_new <- comp_mu_loglik(param = param, y=y, xx= X, 
                           offset= offset, summax)
  ll_old <- min(ll_new/2, ll_new*2)
  iter <- 0
  while (abs((ll_old-ll_new)/ll_new)>tol && iter <= 100){
    iter <- iter+1
    ll_old <- ll_new
    betaold <- param[1:q1]
    etaold <- t(X%*%betaold)[1,]
    muold <-  exp(etaold+offset)
    lambdaold <- param[(q1+1):(q1+n)]
    nuold <- param[(q1+n+1):(q1+2*n)]
    gammaold <-  param[(q1+2*n+1):(q1+2*n+q2)]
    log.Z <- logZ(log(lambdaold), nuold, summax = summax)
    ylogfactorialy <- comp_mean_ylogfactorialy(lambdaold, nuold, log.Z, summax)
    logfactorialy <- comp_mean_logfactorialy(lambdaold, nuold, log.Z, summax)
    variances <- comp_variances(lambdaold, nuold, log.Z, summax)
    variances_logfactorialy <- 
      comp_variances_logfactorialy(lambdaold, nuold, log.Z, summax)
    W <- diag(muold^2/variances)
    z <- etaold + (y-muold)/muold
    beta <- solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%z
    Aterm <- (ylogfactorialy- muold*logfactorialy)
    update_score <- ((Aterm*(y-muold)/variances-lgamma(y+1)+logfactorialy)*nuold)%*%S
    update_info_matrix <- matrix(0, nrow=q2, ncol=q2)
    for (i in 1:n){
      update_info_matrix <- update_info_matrix + 
        ((-(Aterm[i])^2/variances[i] +variances_logfactorialy[i])*nuold[i]^2)*S[i,]%*%t(S[i,])
    }
    gamma <- gammaold + update_score%*%solve(update_info_matrix)
    nu <- exp(t(S%*%t(gamma))[1,])
    eta <- t(X%*%beta)[1,]
    mu <- exp(eta+offset)
    lambda.ok <- comp_lambdas(mu, nu, lambdalb = lambdalb, 
                              lambdaub = min(lambdaub, max(lambdaold*2)), 
                              maxlambdaiter = maxlambdaiter, tol = tol,
                              lambdaint = lambdaold, summax = summax)
    lambda <- lambda.ok$lambda
    lambdaub <-lambda.ok$lambdaub
    param <- c(beta, lambda, nu, gamma) 
    ll_new <- comp_mu_loglik(param = param, y=y, xx= X, 
                             offset= offset, summax)
    nhalf = 0 
    while (ll_new < ll_old && nhalf <=20 && max(gamma-gammaold)>1e-6){
      nhalf <- nhalf+1
      beta <- (beta+betaold)/2
      eta <- t(X%*%beta)[1,]
      mu <- exp(eta+offset)
      gamma <- (gamma + gammaold)/2
      nu <- exp(t(S%*%t(gamma))[1,])
      lambda.ok <- comp_lambdas(mu, nu, lambdalb = lambdalb, lambdaub = lambdaub, 
                                maxlambdaiter = maxlambdaiter, tol = tol,
                                lambdaint = lambdaold, summax = summax)
      lambda <- lambda.ok$lambda
      lambdaub <-lambda.ok$lambdaub
      param <- c(beta, lambda, nu, gamma) 
      ll_new <- comp_mu_loglik(param = param, y=y, xx= X, 
                               offset= offset, summax)
    }
  }
  maxl <-  ll_new # maximum loglikelihood achieved
  beta <-  param[1:q1]  # estimated regression coefficients beta
  lambda <-  param[(q1+1):(q1+n)]  # estimated rates (not generally useful)
  nu <-  param[(q1+n+1):(q1+2*n)] # estimate dispersion
  gamma <-  param[(q1+2*n+1):(q1+2*n+q2)]
  log.Z <- logZ(log(lambda), nu, summax = summax)
  variances <- comp_variances(lambda, nu, log.Z = log.Z, summax = summax)
  eta <- t(X%*%beta)[1,]
  if (is.null(offset)){
    fitted <- exp(eta)
  } else {
    fitted <- exp(eta+offset)
  }
  precision_beta <-  0
  for (i in 1:n){
    precision_beta <-  precision_beta + fitted[i]^2*X[i,]%*%t(X[i,])/variances[i]
  }
  variance_beta <- solve(precision_beta)
  se_beta <- as.vector(sqrt(diag(variance_beta)))
  ylogfactorialy <- comp_mean_ylogfactorialy(lambda, nu, log.Z, summax)
  logfactorialy <- comp_mean_logfactorialy(lambda, nu, log.Z, summax)
  Aterm <- (ylogfactorialy- mu*logfactorialy)
  variances_logfactorialy <- 
    comp_variances_logfactorialy(lambda, nu, log.Z, summax)
  update_info_matrix <- matrix(0, nrow=q2, ncol=q2)
  for (i in 1:n){
    update_info_matrix <- update_info_matrix + 
      ((-(Aterm[i])^2/variances[i] + variances_logfactorialy[i])*nu[i]^2)*(S[i,]%*%t(S[i,]))
  }
  variance_gamma <- solve(update_info_matrix)
  se_gamma <- as.vector(sqrt(diag(variance_gamma)))
  Xtilde <- diag(fitted/sqrt(variances))%*%as.matrix(X)
  h <- diag(Xtilde%*%solve(t(Xtilde)%*%Xtilde)%*%t(Xtilde))
  df.residuals <- length(y)-length(beta)-length(gamma)
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
  out$const_nu <- FALSE
  out$y <- y
  out$x <- X
  out$nobs <- n
  out$iter <- iter
  out$family <- 
    structure(list(family = "CMP(mu, nu)", link = 'log'), 
              class = "family")
  out$coefficients <- c(beta, gamma)
  out$coefficients_beta <- beta
  out$rank <- length(beta)
  out$lambda <- lambda
  out$log_Z <- log.Z
  out$summax <- summax 
  out$offset <- offset
  out$nu <- nu
  out$coefficients_gamma <- gamma
  out$rank_nu <- length(gamma)
  out$s <- S
  out$lambdaub <- lambdaub
  out$linear_predictors <- eta
  out$maxl <- maxl
  out$fitted_values <- fitted
  out$residuals <- y - fitted
  out$leverage <- h
  names(out$leverage) <- 1:n
  out$d_res <- d.res
  out$variance_beta <- variance_beta
  colnames(out$variance_beta) <- row.names(variance_beta)
  out$variance_gamma <- variance_gamma
  colnames(out$variance_gamma) <- row.names(variance_gamma)
  out$se_beta <- se_beta
  out$se_gamma <- se_gamma
  out$df_residuals <- df.residuals
  out$df_null <- n-1
  out$null_deviance <- 2*(sum(indsat.deviance) -
                            sum(dcomp(y, mu = mean(y), nu = nu, log.p = TRUE, 
                                      summax=summax, lambdalb = min(lambdalb),
                                      lambdaub = max(lambdaub))))
  out$deviance <- out$residual_deviance <- 
    2*(sum(indsat.deviance) -
         sum(dcomp(y, lambda = lambda, nu = nu, 
                   log.p = TRUE, summax=summax)))
  names(out$coefficients_beta) <-  labels(X)[[2]]
  names(out$coefficients_gamma) <-  labels(S)[[2]]
  names(out$coefficients) <- c(paste0("beta_", names(out$coefficients_beta)), paste0("gamma_", names(out$coefficients_gamma)))
  class(out) <- "cmp"
  return(out)
}