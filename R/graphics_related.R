#' PIT Plots for a CMP Object
#' 
#' Two plots for the non-randomized PIT are currently available for checking the 
#' distributional assumption of the fitted CMP model: the PIT histogram, and 
#' the uniform Q-Q plot for PIT.
#' 
#' @usage 
#' histcompPIT(object, bins = 10, line = TRUE, colLine = "red", colHist = "royal blue", 
#' lwdLine = 2, main = NULL, ...)
#' qqcompPIT(object, bins = 10, col1 = "red", col2 = "black", lty1 = 1, 
#' lty2 = 2, type = "l", main = NULL, ...)
#' 
#' @param object an object class "cmp", obtained from a call to \code{glm.cmp}.
#' @param bins numeric: the number of bins shown in the PIT histogram or the 
#' PIT Q-Q plot. 
#' @param line logical: if \code{TRUE} (default), the line for displaying the standard 
#' uniform distribution will be shown for the purpose of comparison. 
#' @param colLine	 numeric or charater: the colour of the line for comparison 
#' in PIT histogram.
#' @param lwdLine numeric: the line widths for the comparison line in PIT histogram.
#' @param colHist numeric or character: the colour of the histogram for PIT.
#' @param col1 numeric or character: the colour of the sample uniform Q-Q plot in PIT.
#' @param col2 numeric or character: the colour of the theoretical uniform Q-Q plot in PIT.
#' @param lty1 integer or character string: the line types for the sample 
#' uniform Q-Q plot in PIT, see par(lty = .).
#' @param lty2 an integer or character string: the line types for the theoretical uniform 
#' Q-Q plot in PIT, see par(lty = .).
#' @param type 1-character string: the type of plot for the sample uniform Q-Q plot in PIT.
#' @param main character string: a main title for the plot. 
#' @param ... other arguments passed to plot.default and plot.ts.
#' @details 
#' The histogram and the Q-Q plot are used to compare the fitted profile with a standard
#' uniform distribution. If they match relatively well, it means the CMP distribution 
#' is appropriate for the data. 
#' @references 
#' Czado, C., Gneiting, T. and Held, L. (2009) Predictive model assessment
#' for count data. \emph{Biometrics}, \strong{65}, 1254--1261.
#' @examples 
#' For examples see example(plot.cmp)

histcompPIT <- function (object, bins = 10, line = TRUE, colLine = "red", colHist = "royal blue", lwdLine = 2, main = NULL, ...)
{
  PIT <- compPIT(object, bins = bins)$PIT
  height <- diff(PIT[, ncol(PIT)]) * bins
  if (max(height) > 2) {
    y.upper <- max(height) + 1/(bins/2)
  }
  else {
    y.upper <- 2
  }
  if (is.null(main)) {
    main <- paste("PIT for COM-Poisson")
  }
  barplot(height, ylim = c(0, y.upper), border = TRUE, space = 0,
          xpd = FALSE, xlab = "Probability Integral Transform",
          ylab = "Relative Frequency", main = main, col = colHist,
          ...)
  if (line == TRUE) {
    abline(h = 1, lty = 2, col = colLine, lwd = lwdLine)
  }
  plot.window(xlim = c(0, 1), ylim = c(0, y.upper))
  axis(1)
  box()
}

#' PIT Plots for a CMP Object
#' 
#' Two plots for the non-randomized PIT are currently available for checking the 
#' distributional assumption of the fitted CMP model: the PIT histogram, and 
#' the uniform Q-Q plot for PIT.
#' 
#' @usage 
#' histcompPIT(object, bins = 10, line = TRUE, colLine = "red", colHist = "royal blue", 
#' lwdLine = 2, main = NULL, ...)
#' qqcompPIT(object, bins = 10, col1 = "red", col2 = "black", lty1 = 1, 
#' lty2 = 2, type = "l", main = NULL, ...)
#' 
#' @param object an object class "cmp", obtained from a call to \code{glm.cmp}.
#' @param bins numeric: the number of bins shown in the PIT histogram or the 
#' PIT Q-Q plot. 
#' @param line logical: if \code{TRUE} (default), the line for displaying the standard 
#' uniform distribution will be shown for the purpose of comparison. 
#' @param colLine	 numeric or charater: the colour of the line for comparison 
#' in PIT histogram.
#' @param lwdLine numeric: the line widths for the comparison line in PIT histogram.
#' @param colHist numeric or character: the colour of the histogram for PIT.
#' @param col1 numeric or character: the colour of the sample uniform Q-Q plot in PIT.
#' @param col2 numeric or character: the colour of the theoretical uniform Q-Q plot in PIT.
#' @param lty1 integer or character string: the line types for the sample 
#' uniform Q-Q plot in PIT, see par(lty = .).
#' @param lty2 an integer or character string: the line types for the theoretical uniform 
#' Q-Q plot in PIT, see par(lty = .).
#' @param type 1-character string: the type of plot for the sample uniform Q-Q plot in PIT.
#' @param main character string: a main title for the plot. 
#' @param ... other arguments passed to plot.default and plot.ts.
#' @details 
#' The histogram and the Q-Q plot are used to compare the fitted profile with a standard
#' uniform distribution. If they match relatively well, it means the CMP distribution 
#' is appropriate for the data. 
#' @references 
#' Czado, C., Gneiting, T. and Held, L. (2009) Predictive model assessment
#' for count data. \emph{Biometrics}, \strong{65}, 1254--1261.
#' Jung, R.C and Tremayne, A.R (2011) Useful models for time series of counts or simply wrong ones? \emph{Advances in Statistical Analysis}, \strong{95}, 59–-91.
#' @examples 
#' For examples see example(plot.cmp)

qqcompPIT <- function(object, bins = 10, col1 = "red", col2 = "black", lty1 = 1,
                      lty2 = 2, type = "l", main = NULL, ...){
  dummy.variable <- seq(0, 1, by = 1/bins)
  PIT <- compPIT(object, bins = bins)$PIT
  qq.plot <- PIT[, ncol(PIT)]
  if (is.null(main)) {
    main <- "Q-Q Plot of Uniform PIT"}
  plot(dummy.variable, qq.plot, lty = lty1, col = col1, xlim = c(0, 1), ylim = c(0, 1), type = type, xlab = "Theoretical", ylab = "Sample", main = main, ...)
  abline(0, 1, col = col2, lty = lty2)
  list(sample = qq.plot, theoretical = dummy.variable)
  invisible()
}

#' Non-randomized Probability Integral Transform
#' 
#' Functions to produce the non-randomized probability integral transform (PIT) to check the
#' adequacy of the CMP distribution. 
#' @usage 
#' compPIT(object, bins = 10)
#' @param object an object class "cmp", obtained from a call to \code{glm.cmp}.
#' @param bins numeric: the number of bins used in the PIT. 
#' @details 
#' These functions are used for the assessment of predictive distributions 
#' in discrete data. They obtain the predictive probabilities and the probability 
#' integral transformation for a fitted CMP model.
#' @return
#' compPredProb returns a list with values:
#' \item{upper}{the predictive cumulative probabilities used as the upper bound for computing the non-randomized PIT.}
#' \item{lower}{the predictive cumulative probabilities used as the lower bound for computing the non-randomized PIT.}
#' 
#' compPIT returns a list with values:
#' \item{upper}{the predictive cumulative probabilities used as the upper bound for 
#' computing the non-randomized PIT.}
#' \item{lower}{the predictive cumulative probabilities used as the lower bound for 
#' computing the non-randomized PIT.}
#' \item{conditionalPIT}{the conditional probability integral transformation 
#' given the observed counts.}
#' \item{PIT}{the probability integral transformation.}
#' @references 
#' Czado, C., Gneiting, T. and Held, L. (2009) Predictive model assessment
#' for count data. \emph{Biometrics}, \strong{65}, 1254--1261.
#' Jung, R.C and Tremayne, A.R (2011) Useful models for time series of counts or simply wrong ones? \emph{Advances in Statistical Analysis}, \strong{95}, 59–-91.
#' @examples 
#' 
compPIT <- function (object, bins = 10)
{
  dummy.variable <- seq(0, 1, by = 1/bins)
  predProb <- compPredProb(object)
  con.PIT <- matrix(0, ncol = length(object$y), nrow = length(dummy.variable))
  for (i in 1:length(object$y)) {
    for (j in 1:length(dummy.variable)) {
      if (dummy.variable[j] <= predProb$lower[i])
        con.PIT[j, i] = 0
      if (dummy.variable[j] > predProb$lower[i] & dummy.variable[j] <
          predProb$upper[i]) {
        con.PIT[j, i] = (dummy.variable[j] - predProb$lower[i])/(predProb$upper[i] -
                                                                   predProb$lower[i])
      }
      if (dummy.variable[j] >= predProb$upper[i])
        con.PIT[j, i] = 1
    }
  }
  PIT <- matrix(0, ncol = length(object$y), nrow = length(dummy.variable))
  for (i in 2:length(object$y)) {
    PIT[, i] <- apply(con.PIT[, 1:i], 1, sum) * (i - 1)^(-1)
  }
  list(ConditionalPIT = con.PIT, PIT = PIT)
}

compPredProb <- function (object) {
  lower <- pcomp(object$y-1, nu = object$nu, lambda = object$lambda)
  upper <- lower + dcomp(object$y, nu = object$nu, lambda = object$lambda)
  list(lower = lower, upper = upper)
}


compnormRandPIT <- function (object) {
  temp <- compPredProb(object)
  rt <- qnorm(runif(length(temp$lower), temp$lower, temp$upper))
  rtMid <- qnorm((temp$lower + temp$upper)/2)
  list(rt = rt, rtMid = rtMid)
}


plot.cmp <- function(object, which=c(1L,2L,6L,8L), ask = prod(par("mfcol")<length(which)) && dev.interactive(), bins=10){
  # plot 1 deviance residuals vs fitted
  # plot 2 Histogram of non-randomized PIT
  # plot 3 q-q plot of non-randomized PIT
  # plot 4 Histogram of Randomized Residuals
  # plot 5 Q-Q Plot of Randomized Residuals
  # plot 6 sclae-location plot
  # plot 7 cook's distance vs obs number
  # plot 8 std pearson resid. vs leverage
  show <- rep(FALSE, 9)
  show[which] <- TRUE
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  if (any(show[c(1L,6L)] == TRUE)) {
    dev.res <- object$d.res
  }
  if (show[1L]){
    dev.hold()
    plot(object$linear.predictors, dev.res,
         main= "Residuals vs Fitted",
         ylab = "(Deviance) Residuals",
         xlab = paste(c("Linear Predicted values", object$call)))
    abline(h=0,lty=2)
    lines(lowess(object$linear.predictors,dev.res), col="red")
    index.dev.res <- order(abs(dev.res), decreasing = TRUE)[1:3]
    text(object$linear.predictors[index.dev.res],
         dev.res[index.dev.res], labels=paste(index.dev.res), cex = 0.7, pos = 4)
    dev.flush()
  }
  if (show[2L]){
    dev.hold()
    histcompPIT(object)
    dev.flush()
  }
  if (show[3L]){
    dev.hold()
    qqcompPIT(object)
    dev.flush()
  }
  if (any(show[4L:5L] == TRUE)) {
    rt <- compnormRandPIT(object)$rt
  }
  if (show[4L]) {
    dev.hold()
    hist(rt, main = "Histogram of Randomized Residuals", col = "royal blue",
         xlab = expression(r[t]))
    box()
    dev.flush()
  }
  if (show[5L]) {
    dev.hold()
    qqnorm(rt, main = "Q-Q Plot of Randomized Residuals")
    abline(0, 1, lty = 2, col = "black")
    dev.flush()
  }
  if (show[6L]){
    dev.hold()
    std.dev.res <- dev.res/sqrt(1-object$leverage)
    res <- sqrt(abs(std.dev.res))
    plot(object$linear.predictors, res,
         main="Scale-Location", xlab = paste(c("Linear Predicted values", object$call)),
         ylab = expression(sqrt("|Std. deviance resid.|")))
    lines(lowess(object$linear.predictors,res), col="red")
    index.dev.res <- order(res, decreasing = TRUE)[1:3]
    text(object$linear.predictors[index.dev.res],
         res[index.dev.res], labels=paste(index.dev.res), cex = 0.7, pos = 4)
    dev.flush()
  }
  if (any(show[7L:8L] == TRUE)) {
    rk <- dim(object$x)[2]
    h <- object$leverage
    pear <- residuals(object,type="pearson")
    std.pear <- pear/sqrt(1 - h)
    cook <- (h * std.pear^2)/((1 - h) * rk)
    n <- length(cook)
    index.cook <- order(cook,decreasing = TRUE)[1:3]
  }
  if (show[7L]) {
    dev.hold()
    plot(cook, type = "h", main = "Cook's distance",
         xlab= paste(c("Obs. number", object$call)),
         ylab = "Approx. Cook's distance", ylim = c(0, ymx <- max(cook) * 1.075))
    text(index.cook, cook[index.cook], labels=paste(index.cook), cex = 0.7, pos = 3)
    dev.hold()
  }
  if (show[8L]){
    ylim <- extendrange(r = range(std.pear, na.rm = TRUE), f = 0.08)
    xlim <- extendrange(r = c(0,max(h, na.rm = TRUE)), f = 0.08)
    plot(object$leverage, std.pear, main="Residuals vs Leverage",
         xlab= paste(c("Leverage", object$call)), ylab ="Std. Pearson resid.",
         xlim = xlim, ylim = ylim)
    abline(h = 0, v = 0, lty = 3, col = "gray")
    bound <- par("usr")
    curve(sqrt(rk*0.5*(1-x)/x), 0.01, bound[2], add = TRUE, col = "red",lty=2)
    curve(sqrt(rk*(1-x)/x), 0.01, bound[2], add = TRUE, col = "red",lty=2)
    curve(-sqrt(rk*0.5*(1-x)/x), 0.01, bound[2], add = TRUE, col = "red",lty=2)
    curve(-sqrt(rk*(1-x)/x), 0.01, bound[2], add = TRUE, col = "red",lty=2)
    legend("bottomleft", legend = "Cook's distance",
           lty = 2, col = 2, bty = "n")
    xmax <- min(0.99, bound[2L])
    ymult <- sqrt(rk * (1 - xmax)/xmax)
    cook.levels <- c(0.5,1)
    aty <- sqrt(cook.levels) * ymult
    axis(4, at = c(-rev(aty), aty), labels = paste(c(rev(cook.levels), cook.levels)),
         mgp = c(0.25, 0.25, 0), las = 2,
         tck = 0, cex.axis = 0.75, col.axis = 2)
    text(h[index.cook], std.pear[index.cook], labels=paste(index.cook), cex = 0.7, pos = 4)
    dev.flush()
  }
  invisible()
}

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

# can use generic update
update <- function (object, ...) {
  UseMethod("update")}

predict.cmp <- function(object, newdata = NULL, se.fit = FALSE, type = c("link", "response")){
  type <- match.arg(type)
  if (is.null(newdata)){
    pred <- switch(type, link = object$linear.predictors,
                   response = object$fitted.values)
    if (se.fit){
      se <- switch(type, link = sqrt(diag(object$x%*%object$variance_beta%*%t(object$x))),
                   response = sqrt(diag(object$x%*%object$variance_beta%*%t(object$x)))*
                     object$fitted.values)
      pred <- list(fit = pred, se.fit = se)
    }
  } else {
    mf <- model.frame(delete.response(terms(object)), data=newdata)
    X <- model.matrix(delete.response(terms(object)), mf)
    pred <- switch(type, link = X%*%object$coefficients,
                   response = exp(X%*%object$coefficients))
    if (se.fit){
      se <- switch(type, link = sqrt(diag(X%*%object$variance_beta%*%X)),
                   response = sqrt(diag(X%*%object$variance_beta%*%t(X)))*pred)
      pred <- list(fit = pred, se.fit = se)
    }
  }
  return(pred)
}