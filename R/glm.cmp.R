#' Fit a Mean Parametrized Conway-Maxwell Poisson Generalized Linear Model
#'
#' The function \code{glm.cmp} is used to fit a mean parametrized Conway-Maxwell Poisson
#' generalized linear model with a log-link by using Fisher Scoring iteration.
#'
#' @param formula an object of class 'formula': a symbolic description of the model to be fitted to the mean via log-link.
#' @param formula_nu an optional object of class 'formula': a symbolic description of the model to be fitted to the dispersion via log-link.
#' @param data an optional data frame containing the variables in the model
#' @param offset this can be used to specify an *a priori* known component to be included
#' in the linear predictor for mean during fitting. This should be \code{NULL} or a numeric vector
#' of length equal to the number of cases.
#' @param subset an optional vector specifying a subset of observations to be used in the
#' fitting process.
#' @param na.action a function which indicates what should happen when the data contain
#' NAs. The default is set by the na.action setting of options, and is na.fail if that
#' is unset. The ‘factory-fresh’ default is na.omit. Another possible value is NULL,
#' no action. Value na.exclude can be useful.
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
#' @param contrasts_mu,contrasts_nu optional lists. See the contrasts.arg of model.matrix.default.
#' @export
#' @import stats
#' @details
#' Fit a mean-parametrized COM-Poisson regression using maximum likelihood estimation
#' via an iterative Fisher Scoring algorithm.
#'
#' Currently, the COM-Poisson regression model allows constant dispersion and regression being linked to the dispersion parameter i.e. varying dispersion.
#'
#' For the constant dispersion model, the model is
#'
#' \deqn{Y_i ~ CMP(\mu_i, \nu),}
#'
#' where
#'
#' \deqn{E(Y_i) = \mu_i = exp(x_i^T \beta),}
#'
#' and \eqn{\nu > 0} is the dispersion parameter.
#'
#' The fitted COM-Poisson distribution is over- or under-dispersed
#' if \eqn{\nu < 1} and \eqn{\nu > 1} respectively.
#'
#' For the varying dispersion model, the model is
#'
#' \deqn{Y_i ~ CMP(\mu_i, \nu_i),}
#'
#' where
#'
#' \deqn{E(Y_i) = \mu_i = exp(x_i^T \beta),}
#'
#' and dispersion parameters are model via
#'
#' \deqn{\nu_i = exp(s_i^T \gamma),}
#'
#' where \eqn{x_i} and \eqn{s_i} are some covariates.
#' @return
#' A fitted model object of class \code{cmp} similar to one obtained from \code{glm}
#' or \code{glm.nb}.
#'
#' The function \code{summary} (i.e., \code{\link{summary.cmp}}) can be used to obtain
#' and print a summary of the results.
#'
#' The functions \code{plot} (i.e., \code{\link{plot.cmp}}) and
#' \code{gg_plot} can be used to produce a range
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
#' The functions \code{LRTnu} and \code{cmplrtest} can be used to perform a likelihood ratio
#' chi-squared test for nu = 1 and for nested COM-Poisson model respectively. 
#'
#' An object class 'glm.cmp' is a list containing at least the following components:
#'
#' \item{coefficients}{a named vector of coefficients}
#' \item{coefficients_beta}{a named vector of mean coefficients}
#' \item{coefficients_gamma}{a named vector of dispersion coefficients}
#' \item{se_beta}{approximate standard errors (using observed rather than expected information) for mean coefficients}
#' \item{se_gamma}{approximate standard errors (using observed rather than expected information) for dispersion coefficients}
#' \item{residuals}{the \emph{response} residuals (i.e., observed-fitted)}
#' \item{fitted_values}{the fitted mean values}
#' \item{rank_mu}{the numeric rank of the fitted linear model for mean}
#' \item{rank_nu}{the numeric rank of the fitted linear model for dispersion}
#' \item{linear_predictors}{the linear fit for mean on log scale}
#' \item{df_residuals}{the residuals degrees of freedom}
#' \item{df_null}{the residual degrees of freedom for the null model}
#' \item{null_deviance}{The deviance for the null model.
#' The null model will include only the intercept.}
#' \item{deviance; residual_deviance}{The residual deviance of the model}
#' \item{y}{the \code{y} vector used.}
#' \item{x}{the model matrix for mean}
#' \item{s}{the model matrix for dispersion}
#' \item{model_mu}{the model frame for mu}
#' \item{model_nu}{the model frame for nu}
#' \item{call}{the matched call}
#' \item{formula}{the formula supplied for mean}
#' \item{formula_nu}{the formula supplied for dispersion}
#' \item{terms_mu}{the \code{terms} object used for mean}
#' \item{terms_nu}{the \code{terms} object used for dispersion}
#' \item{data}{the \code{data} argument}
#' \item{offset}{the \code{offset} vector used}
#' \item{lambdaub}{the final \code{lambdaub} used}
#'
#' @references
#' Fung, T., Alwan, A., Wishart, J. and Huang, A. (2020). \code{mpcmp}: Mean-parametrized
#' Conway-Maxwell Poisson Regression. R package version 0.3.4.
#'
#' Huang, A. (2017). Mean-parametrized Conway-Maxwell-Poisson regression models for
#' dispersed counts. \emph{Statistical Modelling} \bold{17}, 359--380.
#'
#' @seealso
#' \code{\link{summary.cmp}}, \code{\link{autoplot.cmp}}, \code{\link{plot.cmp}}, \code{\link{fitted.cmp}},
#' \code{\link{residuals.cmp}} and \code{\link{LRTnu}}.
#'
#' Additional examples may be found in \code{\link{fish}},
#'  \code{\link{takeoverbids}}, \code{\link{cottonbolls}}.
#'
#' @examples
#' ### Huang (2017) Page 368--370: Overdispersed Attendance data
#' data(attendance)
#' M.attendance <- glm.cmp(daysabs~ gender+math+prog, data=attendance)
#' M.attendance
#' summary(M.attendance)
#' \donttest{plot(M.attendance) # or autoplot(M.attendance)
#' }
#'
#' ### Ribeiro et al. (2013): Varying dispersion as a function of covariates
#' \donttest{data(sitophilus)
#' M.sit <- glm.cmp(formula = ninsect ~ extract, formula_nu = ~extract, data = sitophilus)
#' summary(M.sit)
#' }
#' 

glm.cmp <- function(formula,
                    formula_nu = NULL,
                    data,
                    offset = NULL,
                    subset,
                    na.action,
                    betastart = NULL,
                    gammastart = NULL,
                    lambdalb = 1e-10,
                    lambdaub = 1000,
                    maxlambdaiter = 1e3,
                    tol = 1e-6,
                    contrasts_mu = NULL,
                    contrasts_nu = NULL) {
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action",
               "offset"),
             names(mf),
             0L)
  if (is.null(formula)) {
    stop("formula for beta must be specified (can not be NULL)")
  }
  if (is.null(formula_nu) & !is.null(gammastart)) {
    stop(
      "formula_nu should be specified (should not be NULL) \n
         given that you have provided the starting values for the estimates"
    )
  }
  if (lambdalb >= lambdaub) {
    stop("lower bound for the search of lambda must be smaller than the upper bound")
  }
  if (missing(data)) {
    data <- environment(formula)
  }
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf_mu <- eval(mf, parent.frame())
  mt_mu <- attr(mf_mu, "terms")
  y <- model.response(mf_mu)
  X <- model.matrix(formula, mf_mu, contrasts_mu)
  if (!is.null(formula_nu)) {
    temp <- formula_nu
    if (length(formula_nu) == 2) {
      formula_nu[3] <- formula_nu[2]
      formula_nu[2] <- formula[2]
    }
    mf_nu <- mf
    mf_nu$formula <- formula_nu
    mf_nu <- eval(mf_nu, parent.frame())
    mt_nu <- attr(mf_nu, "terms")
    S <- model.matrix(mt_nu, mf_nu, contrasts_nu)
    formula_nu <- temp
  } else {
    mt_nu <- mf_nu <- S <- NULL
  }
  offset <- as.vector(model.offset(mf_mu))
  if (is.null(offset)) {
    offset.cmp <-  rep(0, length(y))
  } else {
    offset.cmp <-  model.extract(mf_mu, "offset")
  }
  if (is.null(S)) {
    out <- fit_glm_cmp_const_nu(
      y = y,
      X = X,
      offset = offset.cmp,
      betastart = betastart,
      lambdalb = lambdalb,
      lambdaub = lambdaub,
      maxlambdaiter = maxlambdaiter,
      tol = tol
    )
  } else {
    out <- fit_glm_cmp_vary_nu(
      y = y,
      X = X,
      S = S,
      offset = offset.cmp,
      betastart = betastart,
      gammastart = gammastart,
      lambdalb = lambdalb,
      lambdaub = lambdaub,
      maxlambdaiter = maxlambdaiter,
      tol = tol
    )
  }
  out$call <- call
  out$formula <- formula
  if (is.null(formula_nu)) {
    out$contrasts_mu <-
      out$formula_nu <- out$terms_nu <- out$model_nu <- NA
  } else {
    out$formula_nu <- formula_nu
    out$terms_nu <- mt_nu
    out$model_nu <- mf_nu
    out$contrasts_mu <- attr(S, "contrasts")
  }
  out$data <- data
  out$terms_mu <- mt_mu
  out$model_mu <- mf_mu
  out$contrasts_mu <- attr(X, "contrasts")
  out$na.action <- attr(mf_mu, "na.action")
  return(out)
}