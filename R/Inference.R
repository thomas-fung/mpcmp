#' Likelihood Ratio Test for nested COM-Poisson models
#' 
#' Perform a likelihood ratio chi-squared test between nested COM-Poisson models. 
#' The test statistics is calculated as \emph{2*(llik- llik_0)}. The test statistics 
#' has degrees of freedom \emph{r} where \emph{r} is the difference in the number of 
#' parameters between the full and null models. 
#' 
#' 
#' @param object1 an object class 'cmp', obtained from a call to \code{glm.cmp}
#' @param object2 an object class 'cmp', obtained from a call to \code{glm.cmp}
#' @param digits numeric; minimum number of significant digits to be used for most numbers.
#' @import stats
#' @export
#' @references 
#' Huang, A. (2017). Mean-parametrized Conway-Maxwell-Poisson regression models for 
#' dispersed counts. \emph{Statistical Modelling} \bold{17}, 359--380.
#' @seealso \code{\link{glm.cmp}}, \code{\link{update.cmp}}
#' @examples 
#' 
#' ## Testing for the mean coefficients
#' data(takeoverbids)
#'
#' ## Fit full model 
#' M.bids.full <- glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#'     + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#'     
#' ## Fit null model; without whtknght
#' M.bids.null <- update(M.bids.full, .~.-whtknght)
#'      
#' ## Likelihood ratio test for the nested models
#' cmplrtest(M.bids.full, M.bids.null) # order of objects is not important
#'
#' ## Testing for dispersion coefficients
#' data(sitophilus) 
#'  M.sit.full <- glm.cmp(formula = ninsect ~ extract, formula_nu = ~extract, data = sitophilus)
#'  
#' ## Fit null model; dropping extract from dispersion equation
#'  M.sit.null1 <- update(M.sit.full, formula_nu. = ~1)   
#'  cmplrtest(M.sit.null1, M.sit.full)
#' 
#' ## Fit null model; using constant dispersion specification
#' M.sit.null2 <- update(M.sit.full, formula_nu. = NULL)   
#' cmplrtest(M.sit.null2, M.sit.full)
#' 
cmplrtest = function(object1, object2, digits=3) {
  if (class(object1) != "cmp") {
    stop("object1 must be an S3 object of class cmp.")
  }
  if (class(object2) != "cmp") {
    stop("object2 must be an S3 object of class cmp.")
  }
  if (object1$const_nu != object2$const_nu){
    if (object1$const_nu){
      object1 <- update(object1, formula_nu. = ~ 1)
    } else {
      object2 <- update(object2, formula_nu. = ~ 1)
    }
  }
  if ( !(all(names(object1$coefficients) %in% names(object2$coefficients))) &&
       !all(names(object2$coefficients) %in% names(object1$coefficients))) {
    warning(paste0("Neither models' coefficient names are subset of the other. ",
                   "Please make sure the models are nested."))
  }
  if (object1$nobs != object2$nobs) {
    stop("The models have a different number of observations.")
  }
  L1 <- object1$maxl
  L2 <- object2$maxl
  df <- length(object1$coefficients) - length(object2$coefficients)
  ttest <- 2*(L1 - L2)
  if (df<0){
    ttest <- -ttest
    df <- -df
  }
  pval <- 1-pchisq(ttest,df)
  if (pval < 2e-16) {
    pval <- "< 2e-16"
  } else {
    pval <- signif(pval, digits)
  }
  cat("\nLikelihood ratio test for testing both COM-Poisson models are equivalent\n")
  cat("LRT-statistic: ", signif(ttest, digits), "\n")
  cat("Chi-sq degrees of freedom: ", df, "\n")
  cat("P-value: ", pval, "\n")
}

#' Likelihood Ratio Test for nu = 1 of a COM-Poisson model
#' 
#' Perform a likelihood ratio chi-squared test for nu = 1 of a COM-Poisson model. 
#' The test statistics is calculated as \emph{2*(llik- llik_0)} where \emph{llik} and 
#' \emph{llik_0} are the log-likelihood of a COM-Poisson and Poisson model respectively.  
#' The test statistic has 1 degrees of freedom. 
#' 
#' @param object an object class 'cmp', obtained from a call to \code{glm.cmp}
#' @param digits numeric; minimum number of significant digits to be used for most numbers.
#' 
#' @import stats
#' @export
#' @references 
#' Huang, A. (2017). Mean-parametrized Conway-Maxwell-Poisson regression models for 
#' dispersed counts. \emph{Statistical Modelling} \bold{17}, 359--380.
#' @examples 
#' data(takeoverbids)
#' M.bids <- glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#'     + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#' LRTnu(M.bids)
LRTnu <- function(object, digits = 3){
  L1 <- object$maxl
  L2 <- as.vector(logLik(glm(object$formula, data = object$data,
                             family = poisson())))
  ttest <- 2*(L1 - L2)
  pval <- 1-pchisq(ttest,1)
  if (pval < 2e-16) {
    pval <- "< 2e-16"
  }
  else {
    pval <- signif(pval, digits)
  }
  cat("\nLikelihood ratio test for testing nu=1:\n\n")
  cat("Log-Likelihood for Mean-CMP: ", signif(L1, digits), "\n")
  cat("Log-Likelihood for Poisson: ", signif(L2, digits), "\n")
  cat("LRT-statistic: ", signif(ttest, digits), "\n")
  cat("Chi-sq degrees of freedom: ", 1, "\n")
  cat("P-value: ", pval, "\n")
}


#' Update and Re-fit a COM-Poisson Model
#' 
#' \code{update} (i.e., \code{update.cmp}) will update and (by-default) re-fit a model. It is
#' identical to \code{update} in the \code{stats} package. 
#' 
#' @param object an object class 'cmp', obtained from a call to \code{glm.cmp}.
#' @param formula. changes to the existing formula in \code{object} -- see \code{update.formula}
#' @param formula_nu. changes to the existing formula_nu in \code{object} -- see \code{update.formula} for details. It also accepts NULL to not regressing on the dispersion. 
#' @param ... other arguments passed to or from other methods  (currently unused).
#' @param evaluate logical; if \code{TRUE} evaluate the new call otherwise simply return 
#' the call
#' @import stats
#' @export
#' @seealso \code{\link{glm.cmp}}, \code{\link{update.formula}}, \code{\link{cmplrtest}}.
#' 
#' @examples 
#' 
#'# To update the mean regression formula
#' data(takeoverbids)
#'
#' ## Fit full model 
#' M.bids.full <- glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#'     + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#' M.bids.full
#'         
#' ## Dropping whtknght
#' M.bids.null <- update(M.bids.full, .~.-whtknght)
#' M.bids.null
#' 
#' ## To update the dispersion regression formula
#' data(sitophilus)
#'
#' ## Fit full model 
#' M.sit.full <- glm.cmp(formula = ninsect ~ extract, formula_nu = ~extract, data = sitophilus)
#' M.sit.full
#'         
#' ## Dropping extract from the dispersion regression
#' M.sit.null1 <- update(M.sit.full, formula_nu =  ~.-extract)
#' M.sit.null1
#' 
#' ## To not regress on the dispersion at all
#' M.sit.null2 <- update(M.sit.full, formula_nu = NULL)
#' M.sit.null2

update.cmp <- function(object, formula., formula_nu., ..., 
                       evaluate = TRUE) {
  if (is.null(call <- getCall(object))) 
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) 
    call$formula <- update.formula(formula(object), formula.)
  if (!missing(formula_nu.)){
    if (!is.null(formula_nu.)){
      if (!inherits(object$formula_nu, "formula")){
        call$formula_nu <- update.formula(formula(~1), formula_nu.)
      } else {
      call$formula_nu <- update.formula(formula(object$formula_nu),
                                      formula_nu.)
      }
    } else {
      call$formula_nu <- NULL
    }
  }
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate) 
    eval(call, parent.frame())
  else call
}


#' Confidence Intervals for CMP Model Parameters
#'
#' Computes confidence intervals for one or more parameters in a
#' fitted model.
#' @param object an object class 'cmp', obtained from a call to \code{\link{glm.cmp}}.
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names (comparing to those provided by \code{\link{coef}()}) . If missing, all parameters are considered. 
#' @param level the confidence level required.
#' @param ... other arguments passed to or from other methods  (currently unused).
#'
#' @return
#' A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2 in % (by default 2.5% and 97.5%).
#' @export
#'
#' @examples
#' data(attendance)
#' M.attendance <- glm.cmp(daysabs~ gender+math+prog, data=attendance)
#' confint(M.attendance)
#' confint(M.attendance, parm = "math", level = 0.9)
#' 
confint.cmp <- function(object, parm, level = 0.95, ...){
  cf <- coef(object)
  if (object$const_nu){
    ses <- object$se_beta
  } else {
    ses <- c(object$se_beta, object$se_gamma)
  }
  pnames <- names(ses) <- names(cf)
  if (missing(parm)) 
    parm <- pnames
  else if (is.numeric(parm)) 
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  pct <- format.perc(a, 3)
  ci <- array(NA_real_, dim = c(length(parm), 2L), 
              dimnames = list(parm, pct))
  ci[] <- cf[parm] + ses[parm] %o% fac
  ci
}
