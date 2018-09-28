#' PIT Plots for a CMP Object
#' 
#' Two plots for the non-randomized PIT are currently available for checking the 
#' distributional assumption of the fitted CMP model: the PIT histogram, and 
#' the uniform Q-Q plot for PIT. 
#' 
#' @param object an object class "cmp", obtained from a call to \code{glm.cmp}.
#' @param bins numeric; the number of bins shown in the PIT histogram or the 
#' PIT Q-Q plot. 
#' @param line logical; if \code{TRUE} (default), the line for displaying the standard 
#' uniform distribution will be shown for the purpose of comparison. 
#' @param colLine	 numeric or charater: the colour of the line for comparison 
#' in PIT histogram.
#' @param lwdLine numeric; the line widths for the comparison line in PIT histogram.
#' @param colHist numeric or character; the colour of the histogram for PIT.
#' @param col1 numeric or character; the colour of the sample uniform Q-Q plot in PIT.
#' @param col2 numeric or character; the colour of the theoretical uniform Q-Q plot in PIT.
#' @param lty1 integer or character string: the line types for the sample 
#' uniform Q-Q plot in PIT, see par(lty = .).
#' @param lty2 an integer or character string: the line types for the theoretical uniform 
#' Q-Q plot in PIT, see par(lty = .).
#' @param type 1-character string; the type of plot for the sample uniform Q-Q plot in PIT.
#' @param main character string; a main title for the plot. 
#' @param ... other arguments passed to plot.default and plot.ts.
#' @import graphics
#' @details 
#' The histogram and the Q-Q plot are used to compare the fitted profile with a standard
#' uniform distribution. If they match relatively well, it means the CMP distribution 
#' is appropriate for the data. 
#' @references 
#' Czado, C., Gneiting, T. and Held, L. (2009). Predictive model assessment
#' for count data. \emph{Biometrics}, \strong{65}, 1254--1261.
#' 
#' Dunsmuir, W.T.M. and Scott, D.J. (2015). The \code{glarma} Package for Observation-Driven
#' Time Seires Regression of Counts. \emph{Journal of Statistical Software}, 
#' \strong{67}, 1--36. 
#' @examples 
#' ## For examples see example(plot.cmp)
#' @name PIT_Plot
NULL

#' @rdname PIT_Plot
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

#' @rdname PIT_Plot
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
#' Functions to produce the non-randomzied probability integral transform (PIT) to check the 
#' adequacy of the distributional assumption of the COM-Poisson model. The majority of the 
#' code and descriptions are taken from Dunsmuir and Scott (2015).
#' 
#' @param object an object class "cmp", obtained from a call to \code{glm.cmp}.
#' @param bins numeric; the number of bins shown in the PIT histogram or the 
#' PIT Q-Q plot. 
#' @details 
#' These functions are used to obtain the predictive probabilities and the probability 
#' integral transform for a fitted COM-Poisson model. The majority of the code and 
#' descriptions are taken from Dunsmuir and Scott (2015).
#' @return 
#' \code{compPredprob} returns a list with values: 
#' \item{upper}{the predictive cumulative probabilities used as the upper bound for 
#' computing the non-randomzied PIT.}
#' \item{lower}{the predictive cumulative probabilities used as the upper bound for 
#' computing the non-randomzied PIT.}
#' 
#' \code{compPIT} returns a list with values:
#' \item{conditionalPIT}{the conditional probability integral transoformation given the 
#' observed counts.}
#' \item{PIT}{the probability integral transformation.}
#' @references 
#' Czado, C., Gneiting, T. and Held, L. (2009). Predictive model assessment
#' for count data. \emph{Biometrics}, \strong{65}, 1254--1261.
#' 
#' Dunsmuir, W.T.M. and Scott, D.J. (2015). The \code{glarma} Package for Observation-Driven
#' Time Seires Regression of Counts. \emph{Journal of Statistical Software}, 
#' \strong{67}, 1--36. 
#' @examples 
#' data(takeoverbids)
#' M.bids <- glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#'     + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#' compPredProb(M.bids)
#' compPIT(M.bids)
#' @name nrPIT
NULL

#' @rdname nrPIT
#' @export
compPredProb <- function (object) {
  lower <- pcomp(object$y-1, nu = object$nu, lambda = object$lambda)
  upper <- lower + dcomp(object$y, nu = object$nu, lambda = object$lambda)
  list(lower = lower, upper = upper)
}

#' @rdname nrPIT
#' @export
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

#' Random Normal Probability Integral Transform
#' 
#' A function to create the normal conditional (randomized) quantile residuals. 
#' The majority of the code and descriptions are taken from Dunsmuir and Scott (2015).
#' @param object an object class "cmp", obtained from a call to \code{glm.cmp}.
#' 
#' @import stats
#' @details 
#' The function \code{compPredProb} produces the non-randomized probability integral 
#' transform(PIT). It returns estimates of the cumulative predictive probabilities as 
#' upper and lower bounds of a collection of intervals. If the model is correct, a 
#' histogram drawn using these estimated probabilities should resemble a histogram 
#' obtained from a sample from the uniform distribution. 
#' 
#' This function aims to produce observations which instead resemble a sample from 
#' a normal distribution. Such a sample can then be examined by the usual tools for 
#' checking normality, such as histograms and normal Q-Q plots.
#' 
#' For each of the intervals produced by \code{compPredProb}, a random uniform observation 
#' is generated, which is then converted to a normal observation by applying the inverse
#' standard normal distribution function (using \code{qnorm}). The vector of these values 
#' is returned by the function in the list element \code{rt}. In addition non-random 
#' observations which should appear similar to a sample from a normal distribution 
#' are obtained by applying \code{qnorm} to the mid-points of the predictive distribution
#' intervals. The vector of these values is returned by the function in the list element 
#' \code{rtMid}.

#' @return 
#' A list consisting of two elements:
#' \item{rt}{the normal conditioanl randomized quantile residuals}
#' \item{rdMid}{the midpoints of the predictive probability intervals}
#' @references 
#' Berkowitz, J. (2001). Testing density forecasts, with applications to risk management.
#' \emph{Journal of Business \& Economic Statistics}, \bold{19}, 465--474.
#' 
#' Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals. \emph{Journal of
#' Computational and Graphical Statistics}, \bold{5}, 236--244.
#' 
#' Dunsmuir, W.T.M. and Scott, D.J. (2015). The \code{glarma} Package for Observation-Driven
#' Time Seires Regression of Counts. \emph{Journal of Statistical Software}, 
#' \strong{67}, 1--36. 
#' @examples 
#' data(takeoverbids)
#' M.bids <- glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#'     + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#' compnormRandPIT(M.bids)
#' @name rPIT
NULL

#' @rdname rPIT
#' @export
compnormRandPIT <- function (object) {
  temp <- compPredProb(object)
  rt <- qnorm(runif(length(temp$lower), temp$lower, temp$upper))
  rtMid <- qnorm((temp$lower + temp$upper)/2)
  list(rt = rt, rtMid = rtMid)
}

#' Plot Diagnostic for a \code{glm.cmp} Object
#' 
#' Eight plots (selectable by \code{which}) are currently available: 
#' a plot of deviance residuals against fitted values, 
#' a non-randomized PIT histogram, 
#' a uniform Q-Q plot for non-randomzied PIT, 
#' a histogram of the normal randomized residuals, 
#' a Q-Q plot of the normal randomized residuals, 
#' a Scale-Location plot of sqrt(| residuals |) against fitted values
#' a plot of Cook's distances versus row lobels
#' a plot of pearson residauls against leverage. 
#' By default, four plots (number 1, 2, 6, and 8 from this list of plots) are provided. 
#' 
#' @param x an object class 'cmp' object, obtained from a call to \code{glm.cmp}
#' @param which if a subset of plots is required, specify a subset of the numbers 1:8. 
#' See 'Details' below. 
#' @param ask logical; if \code{TRUE}, the user is asked before each plot. 
#' @param bins numeric; the number of bins shown in the PIT histogram or the 
#' PIT Q-Q plot. 
#' @param ... other arguments passed to or from other methods (currently unused).
#'
#' @import stats
#' @import graphics
#' @import grDevices
#' @export
#' @details 
#' The 'Scale-Location' plot, also called 'Spread-Loation' plot, takes the square root of 
#' the absolute standardized deviance residuals (\emph{sqrt|E|}) in order to diminish 
#' skewness is much less skewed than than \emph{|E|} for Gaussian zero-mean E. 
#' 
#' The 'Scale-Location' plot uses the standardized deviance residuals while the 
#' Residual-Leverage plot uses the standardized pearson residuals. They are given as 
#' \emph{R[i]/sqrt(1-h.ii)} where h_ii are the diagonal entries of the hat matrix.  
#' 
#' The Residuals-Leverage plot shows contours of equal Cook's distance for values of 0.5 
#' and 1. 
#' 
#' There are two plots based on the non-randomized probability integral transformation (PIT) 
#' using \code{\link{compPIT}}. These are a histogram and a uniform Q-Q plot. If the 
#' model assumption is appropriate, these plots should reflect a sample obtained 
#' from a uniform distribution. 
#' 
#' There are also two plots based on the normal randomzied residuals calculated 
#' using \code{\link{compnormRandPIT}}. These are a histogram and a normal Q-Q plot. If the model
#' assumption is appropriate, these plots should reflect a sample obtained from a normal
#' distribuiton. 
#'
#' @seealso 
#' \code{\link{compPIT}}, \code{\link{compnormRandPIT}}, \code{\link{glm.cmp}}
#' @examples 
#' data(takeoverbids)
#' M.bids <- glm.cmp(numbids ~ leglrest + rearest + finrest + whtknght 
#'     + bidprem + insthold + size + sizesq + regulatn, data=takeoverbids)
#' 
#' ## The default plots are shown
#' plot(M.bids)
#' 
#' ## The plots for the non-randomzied PIT 
#' # plot(M.bids, which = c(2,3))
plot.cmp <- function(x, which=c(1L,2L,6L,8L), 
                     ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                     bins=10, 
                     ...){
  # plot 1 deviance residuals vs fitted
  # plot 2 Histogram of non-randomized PIT
  # plot 3 q-q plot of non-randomized PIT
  # plot 4 Histogram of Randomized Residuals
  # plot 5 Q-Q Plot of Randomized Residuals
  # plot 6 sclae-location plot
  # plot 7 cook's distance vs obs number
  # plot 8 std pearson resid. vs leverage
  object <- x
  show <- rep(FALSE, 8)
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
    pear <- residuals.cmp(object,type="pearson")
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
    x <- NULL
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
    text(h[index.cook], std.pear[index.cook], labels=paste(index.cook), 
         cex = 0.7, pos = 4)
    dev.flush()
  }
  invisible()
}