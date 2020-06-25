#' The Conway-Maxwell-Poisson (COM-Poisson) Distribution.
#' 
#' Density, distribution function, quantile function and random generation for the 
#' Conway-Maxwell-Poisson distribution with parameter \code{mu} and \code{nu}
#'       
#' @param x,q vector of quantiles 
#' @param p vector of probabilities 
#' @param n number of observations. If \code{length(n)} > 1, the length is taken to 
#' be the number required.
#' @param lambda an alternative way than mu to parametrized the distribution. 
#' Must be strictly positive
#' @param mu,nu mean and dispersion parameters. Must be strictly positive.
#' @param log.p logical; if \code{TRUE}, probabilities/densities \eqn{p} are returned as 
#' \eqn{log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P(X \le x)}, 
#' otherwise, \eqn{P(X>x)}.
#' @param lambdalb,lambdaub numeric: the lower and upper end points for the interval to be
#' searched for lambda(s). 
#' @param maxlambdaiter numeric: the maximum number of iterations allowed to solve 
#' for lambda(s).
#' @param tol numeric: the convergence threshold. A lambda is said to satisfy the 
#' mean constraint if the absolute difference between the calculated mean and mu 
#' is less than tol.
#' @param summax numeric; maximum number of terms to be considered in the truncated sum.
#' @return \code{dcomp} gives the density, \code{pcomp} gives the distribution function, \code{qcomp} gives the quantile function, and \code{rcomp} generates random deviates. 
#' 
#' Invalid arguments will result in return value \code{NaN}, with a warning.
#' 
#' The length of the results is determined by \code{n} for \code{rcomp}, and is the maximum 
#' of the lengths of the numerical arguments for the other functions.
#' 
#' The numerical arguments other than \code{n} are recycled to the length of the results. 
#' Only the first argument of the logical arguments are used. 
#' @examples 
#' dcomp(0:5, mu = 2, nu = 1.2)
#' pcomp(5, mu=2, nu =1.2)
#' p <- (1:9)/10
#' qcomp(p, mu = 2, nu = 0.8)
#' rcomp(10, mu = 2, nu = 0.7)
#' @name COM_Poisson_Distribution
NULL

#' @rdname COM_Poisson_Distribution
#' @export
dcomp <- function(x, mu, nu = 1, lambda, log.p = FALSE, lambdalb = 1e-10, 
                  lambdaub = 1000, maxlambdaiter = 1e3, tol = 1e-6, summax){
  # compute the pmf/density for COMP distirbution with mean mu and dispersion nu
  # x, mu, nu are recycled to match the length of each other.
  # lambdaub will be scaled down/up  if there is 
  # over-/under-dispersion so that the correct lambda can be found 
  if (missing(mu) && missing(lambda)){
    stop('argument "mu" is missing, with no default')
  }
  if (!missing(mu) && !missing(lambda)) {
    stop("specify 'mu' or 'lambda' but not both")
  }
  not.miss.mu <- !missing(mu)
  if (missing(mu)) {
    mu <- Inf
  }
  if (missing(lambda)) {
    lambda <- Inf 
  }
  df <- CBIND(x = x, mu=mu, nu=nu, lambda = lambda)
  x <- df[,1]
  mu <- df[,2]
  nu <- df[,3]
  lambda <- df[,4]
  warn <- FALSE
  if (not.miss.mu) {
    if (missing(summax)){
      summax <- ceiling(max(c(mu+20*sqrt(mu/nu),100)))
    }
    mu.ok.ind <- which(mu>0)
    mu.err.ind <- which(mu <= 0)
    if (length(mu.err.ind)>0){ lambda[mu.err.ind] <- mu[mu.err.ind]}
    if (length(mu.ok.ind)>0){
      lambda.ok <- comp_lambdas(mu[mu.ok.ind], nu[mu.ok.ind], 
                                lambdalb = lambdalb, 
                                lambdaub = lambdaub,
                                #lambdaub = min(lambdaub,2*max(lambdaold))
                                maxlambdaiter = maxlambdaiter, tol = tol, 
                                summax = summax)
      lambda[mu.ok.ind] <- lambda.ok$lambda
      lambdaub <- lambda.ok$lambdaub
    }
  } else {
    #A <- (8*nu^2+12*nu+3)/(96*nu^2*lambda^(1/nu))
    #B <- (1+6*nu)/(144*nu^3*lambda^(2/nu))
    #D <- 1+(nu-1)*(A+B)
    #mu <- rep(max(lambda^(1/nu)-(nu-1)/(2*nu)+1/D*((nu-1)*(-A/nu+2*B/nu))),
    #length(lambda))
    #mu_error <- which(is.nan(mu)>0 | mu< 0)
    mu  <- comp_means(lambda, nu, summax = 500)
    if (missing(summax)){
      summax <- ceiling(max(c(mu+20*sqrt(mu/nu),100)))
      cat("As you do not specify mu nor summax, summax will be calculated based on\n")
      cat("mu which is calcualted by truncated sum at 500.\n")
      cat("If you believe the mean of the distribution is somewhat close to 500,\n")
      cat("you may want to do some experiment with the comp_means() and\n")
      cat("specify summax instead to improve the accuracy.\n")
    }
  }
  # at a vector of yvalues
  pmf <- rep(0,length(x))
  for (i in 1:length(x)) {
    if ((mu[i] == 0 || lambda[i] == 0) && x[i]==0) {
      pmf[i] = 0 # log(1), 1 as the distribution is degenerated at 0 
    } else if (mu[i]< 0 | lambda[i] <0 | nu[i] <=0) {
      pmf[i] <- NaN
      warn <- TRUE
    } else {
      if (!is.wholenumber(x[i])) {
        warning(paste("non-integer x =", x[i]))
        pmf[i] <- -Inf # log(0)
      } else {
        if (x[i]<0){pmf[i]= -Inf } else{ #log(0)
          # pmf <- log(density)
          pmf[i] <- x[i]*log(lambda[i])-(nu[i]*lfactorial(x[i]))-
            logZ(log(lambda[i]), nu[i], summax)
        }
      }
    }
  }
  if (!log.p){ pmf = exp(pmf)}
  if (warn){warning("NaN(s) produced")}
  return(pmf)
}


#' @rdname COM_Poisson_Distribution
#' @export
pcomp <- function(q, mu, nu = 1, lambda, lower.tail = TRUE, log.p = FALSE,
                  lambdalb = 1e-10, lambdaub = 1000, maxlambdaiter = 1e3, tol = 1e-6,
                  summax){
  # compute the distribution function for COMP distirbution with mean mu and dispersion nu
  # q, mu, nu are recycled to match the length of each other;
  # lambdaub will be scaled down/up if there is 
  # over-/under-dispersion so that the correct lambda can be found 
  if (missing(mu) && missing(lambda)){
    stop('argument "mu" is missing, with no default')
  }
  if (!missing(mu) && !missing(lambda)) {
    stop("specify 'mu' or 'lambda' but not both")
  }
  not.miss.mu <- !missing(mu)
  if (missing(mu)) {
    mu <- Inf
  }
  if (missing(lambda)) {
    lambda <- Inf
  }
  
  df <- CBIND(q = q, mu=mu, nu=nu, lambda = lambda)
  q <- df[,1]
  mu <- df[,2]
  nu <- df[,3]
  lambda <- df[,4]
  cdf <- rep(0,length(q))
  warn <- FALSE
  if (not.miss.mu) {
    if (missing(summax)){
      summax <- ceiling(max(c(mu+20*sqrt(mu/nu),100)))
    }
    mu.ok.ind <- which(mu>0)
    mu.err.ind <- which(mu <= 0)
    if (length(mu.err.ind)>0){ lambda[mu.err.ind] <- mu[mu.err.ind]}
    if (length(mu.ok.ind)>0){
      lambda.ok <- comp_lambdas(mu[mu.ok.ind], nu[mu.ok.ind], 
                                lambdalb = lambdalb, lambdaub = lambdaub, 
                                #lambdaub = min(lambdaub,2*max(lambdaold))
                                maxlambdaiter = maxlambdaiter, tol = tol, 
                                summax = summax)
      lambda[mu.ok.ind] <- lambda.ok$lambda
      lambdaub <- lambda.ok$lambdaub
    }
  } else {
    #A <- (8*nu^2+12*nu+3)/(96*nu^2*lambda^(1/nu))
    #B <- (1+6*nu)/(144*nu^3*lambda^(2/nu))
    #D <- 1+(nu-1)*(A+B)
    #mu <- rep(max(lambda^(1/nu)-(nu-1)/(2*nu)+1/D*((nu-1)*(-A/nu+2*B/nu))),
    #length(lambda))
    #mu_error <- which(is.nan(mu)>0 | mu< 0)
    if (missing(summax)){
      mu  <- comp_means(lambda, nu, summax = 500)
      summax <- ceiling(max(c(mu+20*sqrt(mu/nu),100)))
      cat("As you do not specify mu nor summax, summax will be calculated based on\n")
      cat("mu which is calcualted by truncated sum at 500.\n")
      cat("If you believe the mean of the distribution is somewhat close to 500,\n")
      cat("you may want to do some experiment with the comp_means() and\n")
      cat("specify summax instead to improve the accuracy.\n")
    }
  }
  for (i in 1:length(q)) {
    if ( (mu[i] == 0 | lambda[i] ==0) && q[i]>=0){
      cdf[i] = 1
    }
    else if (mu[i]< 0 | lambda[i] < 0 | nu[i] <=0) {
      cdf[i] <- NaN
      warn <- TRUE
    } else {
      if (q[i] >= 0){
        cdf[i] <- sum(dcomp(0:floor(q[i]), nu = nu[i], lambda = lambda[i],
                            summax = summax))
      }
    }
  }
  if (warn){warning("NaNs produced")}
  if (!lower.tail){ cdf = 1-cdf}
  if (log.p){ cdf = log(cdf)}
  return(cdf)
}

#' @rdname COM_Poisson_Distribution
#' @export
qcomp <- function(p, mu, nu = 1, lambda, lower.tail = TRUE, log.p = FALSE,
                  lambdalb = 1e-10, lambdaub = 1000, maxlambdaiter = 1e3, tol = 1e-6,
                  summax){
  # compute the distribution function for COMP distirbution with mean mu and dispersion nu
  # q, mu, nu are recycled to match the length of each other;
  # lambdaub will be halved/doubled if there is over-/under-dispersion so that 
  # the correct lambda can be found 
  if (missing(mu) && missing(lambda)){
    stop('argument "mu" is missing, with no default')
  }
  if (!missing(mu) && !missing(lambda)) {
    stop("specify 'mu' or 'lambda' but not both")
  }
  not.miss.mu <- !missing(mu)
  if (missing(mu)) {
    mu <- Inf
  }
  if (missing(lambda)) {
    lambda = Inf
  }
  df <- CBIND(p = p, mu=mu, nu=nu, lambda = lambda)
  p <- df[,1]
  mu <- df[,2]
  nu <- df[,3]
  lambda <- df[,4]
  q <- rep(0,length(p))
  warn <- FALSE
  if (not.miss.mu) {
    if (missing(summax)){
      summax <- ceiling(max(c(mu+20*sqrt(mu/nu),100)))
    }
    mu.ok.ind <- which(mu>0)
    mu.err.ind <- which(mu <= 0)
    if (length(mu.err.ind)>0){ lambda[mu.err.ind] <- mu[mu.err.ind]}
    if (length(mu.ok.ind)>0){
      lambda.ok <- comp_lambdas(mu[mu.ok.ind], nu[mu.ok.ind], 
                                lambdalb = lambdalb, lambdaub = lambdaub, 
                                #lambdaub = min(lambdaub,2*max(lambdaold))
                                maxlambdaiter = maxlambdaiter, tol = tol, 
                                summax = summax)
      lambda[mu.ok.ind] <- lambda.ok$lambda
      lambdaub <- lambda.ok$lambdaub
    }
  } else {
    #A <- (8*nu^2+12*nu+3)/(96*nu^2*lambda^(1/nu))
    #B <- (1+6*nu)/(144*nu^3*lambda^(2/nu))
    #D <- 1+(nu-1)*(A+B)
    #mu <- rep(max(lambda^(1/nu)-(nu-1)/(2*nu)+1/D*((nu-1)*(-A/nu+2*B/nu))),
    #length(lambda))
    #mu_error <- which(is.nan(mu)>0 | mu< 0)
    if (missing(summax)){
      mu  <- comp_means(lambda, nu, summax = 500)
      summax <- ceiling(max(c(mu+20*sqrt(mu/nu),100)))
      cat("As you do not specify mu nor summax, summax will be calculated based on\n")
      cat("mu which is calcualted by truncated sum at 500.\n")
      cat("If you believe the mean of the distribution is somewhat close to 500,\n")
      cat("you may want to do some experiment with the comp_means() and\n")
      cat("specify summax instead to improve the accuracy.\n")
    }
  }
  if (!lower.tail){ p <- 1-p}
  if (log.p){ p <- exp(p)}
  for (i in 1:length(p)) {
    if (mu[i] == 0 | lambda[i] == 0) {
      q[i] <- 0
    }
    else if (mu[i]< 0 | lambda[i] <0 | nu[i] <=0 |p[i]<0 | p[i]>1) {
      q[i] <- NaN
      warn <- TRUE
    } else {
      y <- 0
      py <- dcomp(y, nu = nu[i], lambda = lambda[i], summax=summax)
      while (py <= p[i]){
        y = y+1
        py <- py + dcomp(y, nu = nu[i], lambda = lambda[i], summax = summax)
      }
      q[i] = y
    }
  }
  if (warn){ warning("NaNs produced") }
  return(q)
}

#' @rdname COM_Poisson_Distribution
#' @export
rcomp <- function(n, mu, nu = 1, lambda, lambdalb = 1e-10, 
                  lambdaub = 1000, maxlambdaiter = 1e3, tol = 1e-6,
                  summax){
  # generates random deviates of CMP variables with mean mu and dispersion nu
  # test to see at least one of mu and lambda is missing
  # mu, nu, lambda are recycled to give vectors length n
  # lambdaub will be scaled down/up if there is 
  # over-/under-dispersion so that the correct lambda can be found 
  if (length(n)>1){
    n <- length(n)
  }
  if (missing(mu) && missing(lambda)){
    stop('argument "mu" is missing, with no default')
  }
  if (!missing(mu) && !missing(lambda)) {
    stop("specify 'mu' or 'lambda' but not both")
  }
  not.miss.mu <- !missing(mu)
  if (missing(mu)) {
    mu = Inf
  }
  if (missing(lambda)) {
    lambda = Inf
  }
  if (n < max(length(mu), length(nu), length(lambda))){
    stop("unused argument in mu or nu or lambda")
  }
  df <- CBIND(x = rep(0,n), mu=mu, nu=nu, lambda = lambda)
  x <- df[,1]
  mu <- df[,2]
  nu <- df[,3]
  lambda <- df[,4]
  unif <- runif(n)
  warn <- FALSE
  if (not.miss.mu) {
    if (missing(summax)){
      summax <- ceiling(max(c(mu+20*sqrt(mu/nu),100)))
    }
    mu.ok.ind <- which(mu>0)
    mu.err.ind <- which(mu <= 0)
    if (length(mu.err.ind)>0){ lambda[mu.err.ind] <- mu[mu.err.ind]}
    if (length(mu.ok.ind)>0){
      lambda.ok <- comp_lambdas(mu[mu.ok.ind], nu[mu.ok.ind], 
                                lambdalb = lambdalb, lambdaub = lambdaub, 
                                #lambdaub = min(lambdaub,2*max(lambdaold))
                                maxlambdaiter = maxlambdaiter, tol = tol,
                                summax = summax)
      lambda[mu.ok.ind] <- lambda.ok$lambda
      lambdaub <- lambda.ok$lambdaub
    }
  } else {
    #A <- (8*nu^2+12*nu+3)/(96*nu^2*lambda^(1/nu))
    #B <- (1+6*nu)/(144*nu^3*lambda^(2/nu))
    #D <- 1+(nu-1)*(A+B)
    #mu <- rep(max(lambda^(1/nu)-(nu-1)/(2*nu)+1/D*((nu-1)*(-A/nu+2*B/nu))),
    #length(lambda))
    #mu_error <- which(is.nan(mu)>0 | mu< 0)
    if (missing(summax)){
      mu  <- comp_means(lambda, nu, summax = 500)
      summax <- ceiling(max(c(mu+20*sqrt(mu/nu),100)))
      cat("As you do not specify mu nor summax, summax will be calculated based on\n")
      cat("mu which is calcualted by truncated sum at 500.\n")
      cat("If you believe the mean of the distribution is somewhat close to 500,\n")
      cat("you may want to do some experiment with the comp_means() and\n")
      cat("specify summax instead to improve the accuracy.\n")
    }
  }
  for (i in 1:n){
    if (mu[i] ==0 | lambda[i] == 0){
      x[i] = 0
    } else if (mu[i]< 0 | lambda[i] <0 | nu[i] <=0) {
      x[i] <- NA
      warn <- TRUE
    } else {
      y <- 0
      dc <- dcomp(0:summax, nu = nu[i], lambda = lambda[i], summax=summax)
      py <- dc[y+1]
      while (py <= unif[i]){
        y <- y+1
        py <- py + dc[y+1]
      }
      x[i] <- y
    }
  }
  if (warn){ warning("NAs produced") }
  return(x)
}
