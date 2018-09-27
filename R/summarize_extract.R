#' Extract COM-Poisson Model Residuals
#' 
#' \code{residuals} is a generic function which extracts model residuals from objects 
#' returned by the modelling function \code{glm.comp}. \code{resid} is an alias for 
#' \code{residuals} . 
#' 
#' @usage 
#' ## S3 method for class 'cmp' 
#' residuals(object, type = c("deviance","pearson","response"), ...)
#' 
#' @param object an object class 'cmp', obtained from a call to \code{glm.cmp}.
#' @param type the \code{type} of residuals which should be returned. The alternatives are:
#' 'deviance' (default), 'pearson' and 'response'. Can be abbreviated. 
#' @param ... other arguments passed to or from other methods  (currently unused).

#' @return 
#' Residuals extracted from the object \code{object}.
#' 
#' @seealso 
#' \code{\link{coef.cmp}}, \code{\link{fitted.cmp}}, \code{\link{glm.cmp}}
#' @examples 
#' data(takeoverbids)
#' #' M.bids = glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#' + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#' residuals(M.bids)
residuals.cmp <- function(object, type = c("deviance","pearson","response"), ...){
  type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
  nu <- object$nu
  lambda <- object$lambda
  res <- switch(type, deviance = object$d.res,
                pearson = (y-mu)/sqrt(comp_variances(lambda, nu)),
                response = y - mu
  )
  return(res)
}

#' Extract the (Maximized) Log-Likelihood from a COM-Poisson Model Fit
#' 
#' An accessor function used to extract the (maximized) log-likelihood from a 'cmp' object. 
#' @param object an object of class 'cmp' object, obtained from a call to \code{glm.cmp}
#' @param ... other arguments passed to or from other methods  (currently unused).
#' @param x an object of class 'logLik.cmp', obtained from a call to \code{logLik.cmp}.
#' 
#' @seealso 
#' \code{\link{coef.cmp}}, \code{\link{fitted.cmp}}, \code{\link{glm.cmp}}
#' 
#' @examples 
#' data(takeoverbids)
#' M.bids = glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#' + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#' logLik(M.bids)

#' @name logLik.cmp
logLik.cmp <- function(object,...)
{ out <- object$maxl
  attr(out, "df") <- length(object$coefficients)+1
  class(out) <- "logLik.cmp"
  return(out)}

#' @rdname logLik.cmp
print.logLik.cmp <- function(x,...){
  cat("'log Lik. ' ", x, " (df=", attr(x,"df"),")",sep="")
}

#' Extract the Number of Observation from a COM-Poisson Model Fit
#' 
#' An accessor function used to extract the number of observation from a 'cmp' object.
#' 
#' @param object an object class 'cmp' object, obtained from a call to \code{glm.cmp}
#' @param ... other arguments passed to or from other methods  (currently unused).
#' @return 
#' The number of observations extracted from the object \code{object}.
#' @seealso 
#' \code{\link{coef.cmp}}, \code{\link{fitted.cmp}}, \code{\link{glm.cmp}}
#' @examples 
#' data(takeoverbids)
#' #' M.bids = glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#' + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#' nobs(M.bids)
nobs.cmp <- function(object, ...)
{ return(length(object$nobs))}

#' Akaike's Information Criterion
#' 
#' A function calculating Akaike's Information Criterion (AIC) based on the log-likelihood
#' value extracted from \code{\link{logLik.cmp}}, according to the formula 
#' \emph{-2*log-likelihood + k*npar}, where \emph{npar} represents the number of parameters 
#' in the fitted model, and \emph{k=2} for the usual AIC or \emph{k=log(n)} (\emph{n} being 
#' the number of observations) for the so-called BIC (Bayesian Information Criterion).
#' 
#' @param object an object class 'cmp' object, obtained from a call to \code{glm.cmp}
#' @param ... other arguments passed to or from other methods  (currently unused).
#' @param k numeric: the \emph{penalty} per parameter to be used; the default
#' k = 2 is the classical AIC.
#' @details 
#' When comparing models fitted by maximum likelihood to the same data, the smaller the AIC or 
#' BIC, the better the fit. 
#' @return 
#' A numeric value with the corresponding AIC (or BIC, or ..., depends on k).
#' @seealso 
#' \code{\link{logLik.cmp}}, \code{\link{nobs.cmp}}, \code{\link{glm.cmp}}
#' @examples 
#' data(takeoverbids)
#' M.bids = glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#' + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#' AIC(M.bids)
AIC.cmp = function(object, ..., k = 2){
  temp <- logLik(object)
  aic <- -2*as.numeric(temp)+k*attr(temp,"df")
  return(aic)
}

#' Extract Fitted Values from a COM-Poisson Model Fit
#' 
#' An accessor function used to extract the fitted values from a 'cmp' object.
#' \code{fitted.values} is an alias for \code{fitted}.
#' 
#' @param object an object class 'cmp' object, obtained from a call to \code{glm.cmp}
#' @param ... other arguments passed to or from other methods  (currently unused).
#' @return 
#' Fitted values \code{mu} extracted from the object \code{object}.
#' @seealso 
#' \code{\link{coef.cmp}}, \code{\link{residuals.cmp}}, \code{\link{glm.cmp}}.
#' #' @examples 
#' data(takeoverbids)
#' #' M.bids = glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#' + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#' fitted(M.bids)
fitted.cmp = function(object, ...){
  return(object$fitted.values)
}

#' Extract Model Coefficients from a COM-Poisson Model Fit
#' 
#' An accessor function used to extract model coefficients from a 'cmp' object.
#' \code{coefficients} is an alias for \code{coef}.
#' 
#' @param object an object class 'cmp' object, obtained from a call to \code{glm.cmp}
#' @param ... other arguments passed to or from other methods  (currently unused).
#' 
#' @return 
#' Coefficients extracted from the object \code{object}.
#' @seealso 
#' \code{\link{fitted.cmp}}, \code{\link{residuals.cmp}}, \code{\link{glm.cmp}}.
#' 
#' @examples 
#' data(takeoverbids)
#' #' M.bids = glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#' + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#' coef(M.bids)
coef.cmp = function(object, ...){
  return(object$coefficients)
}


#' Summarize COM-Poisson Model Fit
#' 
#' \code{summary} and \code{print} method for class \code{cmp}.  
#' 
#' @param object an object class 'cmp', obtained from a call to \code{glm.cmp}.
#' @param digits numeric; minimum number of significant digits to be used for most numbers.
#' @param ... other arguments passed to or from other methods  (currently unused).

#' @details  
#' \code{summary.glm} tries to be smart about formatting the coefficients, standard errors 
#' and gives 'signifiance starts'. The \code{coefficients} component of the result gives the 
#' estimated coefficients and their estimated standard errors, together with their ratio. This 
#' third column is labelled as \code{Z value} as the dispersion is fixed for this family. A
#' forth column gives the two-tailed p-value corresponding to Z value based on the asymptotic 
#' Normal reference distribuiton. 
#' 
#' @seealso 
#' \code{\link{coef.cmp}}, \code{\link{fitted.cmp}}, \code{\link{glm.cmp}}.
#' @examples 
#' data(takeoverbids)
#' #' M.bids = glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#' + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#' M.bids
#' summary(M.bids)
#' @name summary.cmp
summary.cmp <- function(object, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall: ", paste(deparse(object$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  cat("\nDeviance Residuals:" , "\n")
  if (object$df.residual > 5) {
    residuals.dev = setNames(quantile(object$d.res, na.rm = TRUE),
                             c("Min", "1Q", "Median", "3Q", "Max"))
  }
  xx <- zapsmall(residuals.dev, digits + 1L)
  print.default(xx, digits = digits, na.print = "", print.gap = 2L)
  cat("\nLinear Model Coefficients:\n")
  cmat <- cbind(object$coefficients, object$se_beta)
  cmat <- cbind(round(cmat,4), cmat[, 1]/cmat[, 2])
  cmat <- cbind(cmat, 2*pnorm(-abs(cmat[,3])))
  colnames(cmat) <- c("Estimate", "Std.Err", "Z value", "Pr(>|z|)")
  printCoefmat(cmat, digits = 3, P.values = TRUE)
  cat("\n(Dispersion parameter for Mean-CMP(", signif(object$nu, digits), ") family taken to be 1)\n\n",
      apply(cbind(paste(format(c("Null", "Residual"), justify = "right"), "deviance:"),
                  format(unlist(object[c("null.deviance","residuals.deviance")]),
                         digits = max(5L, digits + 1L)), " on",
                  format(unlist(object[c("df.null", "df.residuals")])), 
                  "degrees of freedom\n"),
            1L, paste, collapse = " "), sep = "")
  cat("\nAIC:", format(AIC(object)), "\n\n")
}

#' @rdname summary.cmp
print.cmp <- function(object,...)
{
  cat("\nCall: ", paste(deparse(object$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  cat("\nLinear Model Coefficients:\n")
  print.default(format(signif(object$coefficients,3)), print.gap = 2,quote = FALSE)
  cat("\nDispersion (nu):", object$nu)
  cat("\nDegrees of Freedom:", object$df.null, "Total (i.e. Null); ",
      object$df.residuals, "Residual")
  cat("\nNull Deviance:", object$null.deviance, "\nResidual Deviance:",
      object$residuals.deviance, "\nAIC:", format(AIC(object)), "\n\n")
}