#' Fit a Mean Parametrized Conway-Maxwell Poisson Generalized Linear
#' Model with constant dispersion.
#'
#' @description
#' This is a workhorse function in which glm.cmp to call upon to fit a
#' mean-parametrized Conway-Maxwell Poisson generalized linear model with
#' constant dispersion.
#'
#' @param y the response y vector.
#' @param X the design matrix for regressing the mean
#' @param offset, this can be used to specify an *a priori* known component to be included
#' in the linear predictor for mean during fitting. This should be \code{NULL} or a numeric vector
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
#' @keywords internal
#' @return
#' A fitted model object of class \code{cmp} similar to one obtained from \code{glm} or \code{glm.nb}.
#'
#' @details
#' If \code{X} is not of full column rank, the aliased (linearly dependent)
#' columns are dropped before fitting -- identified via the same pivoted QR
#' decomposition that \code{\link[stats]{glm.fit}} uses -- and the
#' corresponding coefficients are reported as \code{NA}, mirroring the
#' behaviour of \code{\link[stats]{glm}} for rank-deficient designs.
#'
#' @examples
#' ## For examples see example(glm.cmp)
fit_glm_cmp_const_nu <- function(y = y, X = X, offset = offset,
                                 betastart = betastart,
                                 lambdalb = lambdalb, lambdaub = lambdaub,
                                 maxlambdaiter = maxlambdaiter, tol = tol) {
  X_full <- X
  q_full <- ncol(X_full)
  M0 <- stats::glm(y ~ -1 + X + offset(offset),
    start = betastart, family = stats::poisson()
  )
  offset <- M0$offset
  n <- length(y) # sample size
  # Identify the identifiable (full-rank) columns of X using the pivoted QR
  # decomposition already computed by glm.fit() for M0. Any column not among
  # the first `rank_x` pivoted columns is aliased and is dropped from the
  # fit; its coefficient is reported as NA at the end, as glm() would do.
  rank_x <- M0$rank
  keep_x <- sort(M0$qr$pivot[seq_len(rank_x)])
  if (rank_x < q_full) {
    warning(
      "design matrix 'X' is rank-deficient (rank ", rank_x, " out of ",
      q_full, " column(s)); coefficient(s) for the aliased column(s) ",
      "will be reported as NA, as with glm()."
    )
  }
  X <- X_full[, keep_x, drop = FALSE] # full-rank subset used for fitting
  q <- ncol(X) # number of identifiable covariates for mu
  # starting values for optimization
  beta0 <- stats::coef(M0)[keep_x]
  lambda0 <- mu0 <- M0$fitted.values
  nu_lb <- 1e-10
  summax <- ceiling(max(c(max(y) + 20 * sqrt(var(y)), 100)))
  nu0 <- exp(optimize(
    f = comp_mu_neg_loglik_log_nu_only,
    interval = c(max(-log(1 + 2 * mu0)), 5), mu = mu0, y = y,
    summax = summax
  )$minimum)
  lambda0 <- (mu0 + (nu0 - 1) / (2 * nu0))^(nu0)
  summax <- ceiling(max(c(mu0 + 20 * sqrt(mu0 / nu0), 100)))
  lambda.ok <- comp_lambdas(mu0, nu0,
    lambdalb = lambdalb,
    lambdaub = min(lambdaub, 2 * max(lambda0)),
    maxlambdaiter = maxlambdaiter, tol = tol, summax = summax,
    lambdaint = lambda0
  )
  lambda0 <- lambda.ok$lambda
  lambdaub <- lambda.ok$lambdaub
  param <- c(beta0, lambda0, nu0)
  ll_old <- as.numeric(logLik(M0))
  ll_new <- comp_mu_loglik(param = param, y = y, xx = X, offset = offset, summax = summax)
  iter <- 1
  while (abs((ll_new - ll_old) / ll_new) > tol && iter <= 100) {
    iter <- iter + 1
    ll_old <- ll_new
    paramold <- param
    betaold <- param[1:q]
    etaold <- t(X %*% betaold)[1, ]
    muold <- exp(etaold + offset)
    lambdaold <- param[(q + 1):(q + n)]
    nuold <- param[q + n + 1]
    log.Z <- logZ(log(lambdaold), nuold, summax = summax)
    ylogfactorialy <- comp_mean_ylogfactorialy(lambdaold, nuold, log.Z, summax)
    logfactorialy <- comp_mean_logfactorialy(lambdaold, nuold, log.Z, summax)
    variances <- comp_variances(lambdaold, nuold, log.Z, summax)
    variances_logfactorialy <-
      comp_variances_logfactorialy(lambdaold, nuold, log.Z, summax)
    W <- diag(muold^2 / variances)
    z <- etaold + (y - muold) / muold
    beta <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% z
    eta <- t(X %*% beta)[1, ]
    mu <- exp(eta + offset)
    Aterm <- (ylogfactorialy - mu * logfactorialy)
    update_score <- sum(Aterm * (y - mu) / variances - (lgamma(y + 1) - logfactorialy))
    update_info_matrix <- sum(-Aterm^2 / variances + variances_logfactorialy)
    if (update_info_matrix < 0) {
      update_info_matrix <- sum(variances_logfactorialy)
    }
    nu <- nuold + update_score / update_info_matrix
    while (nu < nu_lb) {
      nu <- (nu + nuold) / 2
    }
    lambda.ok <- comp_lambdas(mu, nu,
      lambdalb = lambdalb,
      lambdaub = min(lambdaub, 2 * max(lambdaold)),
      maxlambdaiter = maxlambdaiter, tol = tol,
      lambdaint = lambdaold, summax = summax
    )
    lambda <- lambda.ok$lambda
    lambdaub <- lambda.ok$lambdaub
    param <- c(beta, lambda, nu)
    ll_new <- comp_mu_loglik(param = param, y = y, xx = X, offset = offset, summax = summax)
    halfstep <- 0
    while (ll_new < ll_old && halfstep <= 20 && abs((ll_new - ll_old) / ll_new) > tol) {
      halfstep <- halfstep + 1
      beta <- (beta + betaold) / 2
      nu <- (nu + nuold) / 2
      eta <- t(X %*% beta)[1, ]
      mu <- exp(eta + offset)
      lambdaold <- lambda
      lambda.ok <- comp_lambdas(mu, nu,
        lambdalb = lambdalb,
        lambdaub = min(lambdaub, 2 * max(lambdaold)),
        maxlambdaiter = maxlambdaiter, tol = tol,
        lambdaint = lambda, summax = summax
      )
      lambda <- lambda.ok$lambda
      lambdaub <- lambda.ok$lambdaub
      param <- c(beta, lambda, nu)
      ll_new <- comp_mu_loglik(param = param, y = y, xx = X, offset = offset, summax)
    }
  }
  maxl <- ll_new # maximum loglikelihood achieved
  beta <- param[1:q] # estimated regression coefficients beta (identifiable subset)
  lambda <- param[(q + 1):(q + n)] # estimated rates (not generally useful)
  nu <- param[q + n + 1] # estimate dispersion
  precision_beta <- 0
  eta <- t(X %*% beta)[1, ]
  if (is.null(offset)) {
    fitted <- exp(eta)
  } else {
    fitted <- exp(eta + offset)
  }
  log.Z <- logZ(log(lambda), nu, summax = summax)
  variances <- comp_variances(lambda, nu, log.Z = log.Z, summax = summax)
  for (i in 1:n) {
    precision_beta <- precision_beta + fitted[i]^2 * X[i, ] %*% t(X[i, ]) / variances[i]
  }
  variance_beta <- solve(precision_beta)
  se_beta <- as.vector(sqrt(diag(variance_beta)))
  Xtilde <- diag(fitted / sqrt(variances)) %*% as.matrix(X)
  h <- diag(Xtilde %*% solve(t(Xtilde) %*% Xtilde) %*% t(Xtilde))
  df_residuals <- length(y) - rank_x
  if (df_residuals > 0) {
    indsat_deviance <- dcomp(y,
      mu = y, nu = nu, log.p = TRUE, lambdalb = min(lambdalb),
      lambdaub = lambdaub, maxlambdaiter = maxlambdaiter,
      tol = tol, summax = summax
    )
    indred_deviance <- 2 * (indsat_deviance - dcomp(y,
      nu = nu, lambda = lambda,
      log.p = TRUE, summax = summax
    ))
    d_res <- sign(y - fitted) * sqrt(abs(indred_deviance))
  } else {
    d_res <- rep(0, length(y))
  }
  # Expand the identifiable coefficients/SEs back out to X_full's original
  # column order, inserting NA for any aliased column -- mirrors the
  # coefficient vector returned by glm() for a rank-deficient fit.
  beta_full <- rep(NA_real_, q_full)
  names(beta_full) <- colnames(X_full)
  beta_full[keep_x] <- as.vector(beta)
  se_beta_full <- rep(NA_real_, q_full)
  names(se_beta_full) <- colnames(X_full)
  se_beta_full[keep_x] <- se_beta
  out <- list()
  out$const_nu <- TRUE
  out$y <- y
  # Store the full (un-reduced) design matrix, matching the convention of
  # model.matrix.lm()/model.matrix.glm() -- callers such as predict.cmp()
  # drop the aliased columns themselves, using the NA pattern in
  # out$coefficients, before doing any matrix algebra with it.
  out$x <- X_full
  out$family <-
    structure(list(
      family =
        paste0(
          "CMP(mu, ",
          signif(
            nu,
            max(3, getOption("digits") - 4)
          ),
          ")"
        ), link = "log"
    ), class = "family")
  out$nobs <- n
  out$iter <- iter
  out$coefficients <- beta_full
  out$rank <- rank_x
  out$lambda <- lambda
  out$log_Z <- log.Z
  out$summax <- summax
  out$offset <- offset
  out$nu <- nu
  out$lambdaub <- lambdaub
  out$linear_predictors <- eta
  out$maxl <- maxl
  out$fitted_values <- fitted
  out$residuals <- y - fitted
  out$leverage <- h
  names(out$leverage) <- 1:n
  out$d_res <- d_res
  out$variance_beta <- variance_beta
  colnames(out$variance_beta) <- row.names(variance_beta)
  out$se_beta <- se_beta_full
  out$df_residuals <- df_residuals
  out$df_null <- n - 1
  out$s <- NA
  out$formula_nu <- NA
  out$gamma <- NA
  out$variance_gamma <- NA
  out$se_gamma <- NA
  out$null_deviance <- 2 * (sum(indsat_deviance) -
    sum(dcomp(y,
      mu = mean(y), nu = nu, log.p = TRUE,
      summax = summax, lambdalb = min(lambdalb),
      lambdaub = max(lambdaub)
    )))
  out$deviance <- out$residual_deviance <- 2 * (sum(indsat_deviance) -
    sum(dcomp(y,
      lambda = lambda, nu = nu,
      log.p = TRUE, summax = summax
    )))
  class(out) <- "cmp"
  return(out)
}
