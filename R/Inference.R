#' Likelihood Ratio Test for nested COM-Poisson models
#' 
#' Perform a likelihood ratio chi-sqaured test between nested COM-Poisson models. 
#' The test statistics is calculated as \emph{2*(llik- llik_0)}. The test statistics 
#' has degrees of freedom \emph{r} where \emph{r} is the difference in the number of 
#' parameters between the full and null models. 
#' 
#' Obviously the comparison between two models will only be valid if they are fitted 
#' to the same data set. 
#' 
#' @param object1 an object class 'cmp', obtained from a call to \code{glm.cmp}
#' @param object2 an object class 'cmp', obtained from a call to \code{glm.cmp}
#' @param digits numeric; minimum number of significant digits to be used for most numbers.
#' @import stats
#' @export
#' @references 
#' Huang, A. (2017). Mean-parametrized Conway–Maxwell–Poisson regression models for 
#' dispersed counts. \emph{Statistical Modelling} \bold{17}, 359--380.
#' @seealso \code{\link{glm.cmp}}, \code{\link{update.cmp}}
#' @examples 
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
cmplrtest = function(object1,object2, digits=3) {
  if (class(object1) != "cmp") {
    stop("object1 must be an S3 object of class cmp.")
  }
  if (class(object2) != "cmp") {
    stop("object2 must be an S3 object of class cmp.")
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
  cat("\nLikelihood ratio test for testing both COP-Poisson models are equivalent\n")
  cat("LRT-statistic: ", signif(ttest, digits), "\n")
  cat("Chi-sq degrees of freedom: ", df, "\n")
  cat("P-value: ", pval, "\n")
}

#' Likelihood Ratio Test for nu = 1 of a COM-Poisson model
#' 
#' Perform a likelihood ratio chi-sqaured test for nu = 1 of a COM-Poisson model. 
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
#' Huang, A. (2017). Mean-parametrized Conway–Maxwell–Poisson regression models for 
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
#' \code{update} (i.e., \code{update.cmp}) will upate and (by-default) re-fit a model. It is
#' identical to \code{update} in the \code{stats} package. 
#' 
#' @param object an object class 'cmp', obtained from a call to \code{glm.cmp}.
#' @param formula. changes to the existing formula in \code{object} -- see \code{update.formula}
#' for details
#' @param ... other arguments passed to or from other methods  (currently unused).
#' @param evaluate logical; if \code{TRUE} evaluate the new call otherwise simply ruturn 
#' the call
#' @import stats
#' @export
#' @seealso \code{\link{glm.cmp}}, \code{\link{update.formula}}, \code{\link{cmplrtest}}.
#' 
#' @examples 
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
update.cmp <- function(object, formula., ..., evaluate = TRUE) {
  if (is.null(call <- getCall(object))) 
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) 
    call$formula <- update.formula(formula(object), formula.)
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

