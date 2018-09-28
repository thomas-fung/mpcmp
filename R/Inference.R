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
  cat("\nLikelihood ratio test for testing both CMP models are equivalent\n")
  cat("LRT-statistic: ", signif(ttest, digits), "\n")
  cat("Chi-sq degrees of freedom: ", df, "\n")
  cat("P-value: ", pval, "\n")
}

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
