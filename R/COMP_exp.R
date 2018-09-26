CBIND <- function(..., deparse.level = 1) {
  dots <- list(...)
  len <- sapply(dots, length)
  dots <- lapply(seq_along(dots),
                 function(i, x, len) rep(x[[i]], length.out = len),
                 x = dots, len = max(len))
  do.call(cbind, c(dots, deparse.level = deparse.level))
}

comp_mean_logfactorialy = function(lambda, nu){
  # approximates mean by truncation of Ylog(Y!) for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(lambda=lambda, nu=nu)
  lambda <- df[,1]
  nu <- df[,2]
  if (length(lambda)>1 && length(nu>1) && length(lambda)!= length(nu)){
    stop("lambda, nu must be scalars or vectors of the same length")}
  summax <- 100
  termlim <- 1e-6
  sum1 <- 0
  for (y in 1:summax){
    term <- lgamma(y)*exp(log(lambda^(y-1)) - nu*lgamma(y))
    if (y > 3) {
      if (max(term/sum1, na.rm = TRUE) < termlim){
        break
      }
    }
    sum1 <- sum1 + term
  }
  mean1 <- sum1/Z(lambda, nu)
  return(mean1)
}

comp_mean_ylogfactorialy <- function(lambda, nu){
  # approximates mean by truncation of Ylog(Y!) for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(lambda=lambda, nu=nu)
  lambda <- df[,1]
  nu <- df[,2]
  summax <- 100
  termlim <- 1e-6
  sum1 <- 0
  for (y in 1:summax) {
    term <- (y-1)*lgamma(y)*exp(log(lambda^(y-1)) - nu*lgamma(y))
    if (y > 3) {
      if ( max(term/sum1,na.rm = TRUE) < termlim){
        break
      }
    }
    sum1 <- sum1 + term
  }
  mean1 <- sum1/Z(lambda, nu)
  return(mean1)
}

comp_means <- function(lambda, nu) {
  # approximates mean by truncation of COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(lambda=lambda, nu=nu)
  lambda <- df[,1]
  nu <- df[,2]
  summax <- 100
  termlim <- 1e-6
  sum1 <- 0
  for (y in 1:summax) {
    term <- (y-1)*exp(log(lambda^(y-1)) - nu*lgamma(y))
    if (y > 3) {
      if (max(term/sum1,na.rm = TRUE) < termlim){
        break
      }
    }
    sum1 <- sum1 + term
    mean1 <- sum1/Z(lambda, nu)
  }
  return(mean1)
}


comp_variances <- function(lambda, nu) {
  # approximates normalizing constant by truncation for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(lambda=lambda, nu=nu)
  lambda <- df[,1]
  nu <- df[,2]
  summax <- 100
  termlim <- 1e-6
  sum2 <- 0
  for (y in 1:summax){
    term <- (y-1)^2*exp(log(lambda^(y-1)) - nu*lgamma(y))
    if (y > 3) {
      if (max(term/sum2,na.rm = TRUE) < termlim) {
        break
      }
    }
    sum2 <- sum2 + term
  }
  var1 <- sum2/Z(lambda,nu) - (comp_means(lambda,nu))^2
  return(var1)
}

comp_variances_logfactorialy <- function(lambda, nu) {
  # approximates normalizing constant by truncation for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(lambda=lambda, nu=nu)
  lambda <- df[,1]
  nu <- df[,2]
  summax <- 100
  termlim <- 1e-6
  sum2 <- 0
  for (y in 1:summax){
    term <- (lgamma(y))^2*exp(log(lambda^(y-1)) - nu*lgamma(y))
    if (y > 3) {
      if (max(term/sum2,na.rm = TRUE) < termlim) {
        break
      }
    }
    sum2 <- sum2 + term
  }
  var1 <- sum2/Z(lambda,nu) - (comp_mean_logfactorialy(lambda,nu))^2
  return(var1)
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

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}


# too slow. replaced.
comp_lambda <- function(mu, nu, type = "exact"){
  if (missing(type)){
    type = "exact"
  }
  #for given scale mean mu and scalar dispersion nu, solve for the scalar rate lambda
  f <- function(x,mu,nu,type){ comp_means(x, nu, type)-mu}
  lambda <- uniroot(f, interval = c(1e-10,1), extendInt = "yes", mu=mu,nu= nu, type=type)$root
  return(lambda)
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

comp_lambdas2 <- function(mu, nu, lambdalb = 1e-10, lambdaub = 1200,
                          maxlambdaiter = 1e3, tol = 1e-6, lambdaint = 1){
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
  tol <- 1e-6
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
  while (length(not.converge.ind)>0 && iter <10000){
    still.above.target.ind = which(mean1[not.converge.ind]
                                   >mu[not.converge.ind])
    still.below.target.ind = which(mean1[not.converge.ind]<mu[not.converge.ind])
    lb[not.converge.ind[still.below.target.ind]] =
      lambda[not.converge.ind[still.below.target.ind]]
    ub[not.converge.ind[still.above.target.ind]] =
      lambda[not.convergeind[still.above.target.ind]]
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
  return(lambda)
}



glm.cmp <- function(formula, data, offset = NULL,
                     lambdalb = 1e-10, lambdaub = 1299, maxlambdaiter = 1e3, tol = 1e-6){
  call <- match.call()
  if (is.null(formula)) {
    stop("formula must be specified (can not be NULL)")
  }
  if (lambdalb>=lambdaub) {
    stop("lower bound for the search of lambda must be smaller than the upper bound")
  }
  if (missing(data))
    data <- environment(formula)
  mf <- model.frame(formula, data=data)
  mt <- attr(mf, "terms")
  y <- model.response(mf)
  X <- model.matrix(formula,mf)
  if (is.null(offset)){
    offset.cmp = rep(0,length(y))
  } else {
    offset.cmp = offset
  }
  M0 <- glm(y~-1+X+offset(offset.cmp), family=poisson())
  offset <- M0$offset
  #y = as.vector(M1$y)
  n <- length(y) # sample size
  #xx = as.matrix(X)
  q <- ncol(X)  # number of covariates
  #starting values for optimization
  #dat = data.frame(y,xx)
  beta0 <- coef(M0)
  nu0 <- 1
  nu_lb <- 1e-10
  lambda0 <- fitted.values(M0)
  param <- c(beta0,lambda0,nu0)
  ll_old <- -comp_mu_loglik_exact(param = param, y=y, xx= X, offset=offset)$objective
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
    ll_new <- -comp_mu_loglik_exact(param = param, y=y, xx= X, offset= offset)$objective
    halfstep <- 0
    while (ll_new < ll_old && halfstep <= 20){
      halfstep <- halfstep + 1
      beta <- (beta+betaold)/2
      nu <- (nu+nuold)/2
      eta <- t(X%*%beta)[1,]
      mu <- exp(eta+offset)
      lambda <- comp_lambdas(mu,nu)
      param <- c(beta, lambda, nu)
      ll_new <- -comp_mu_loglik_exact(param = param, y=y, xx= X, offset= offset)$objective
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
