Z <- function(lambda, nu, log.z = FALSE, summax){
  # approximates normalizing constant for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(lambda=lambda, nu=nu)
  lambda <- df[,1]
  nu <- df[,2]
  termlim <- 1e-6
  # zero vector length of same length as lambda
  sum1 <- term <- rep(0, length(lambda))
  not.converge.ind <- 1:length(lambda)
  for(y in 1:summax) {
    term[not.converge.ind] <- exp((y-1)*log(lambda[not.converge.ind]) - 
                                    nu[not.converge.ind]*lgamma(y))
    sum1[is.infinite(term)] <- Inf
    term[is.infinite(sum1)] <- 0
    if (y > 3){
      not.converge.ind = which(term/sum1 >= termlim)
      if (length(not.converge.ind)==0) {
        break
      }
    }
    sum1 <- sum1+term
  }
  sum1 <- log(sum1)
  infinite.ind <- which(is.infinite(sum1))
  if (length(infinite.ind)!=0){
    nu.sub <- nu[infinite.ind]
    lambda.sub <- lambda[infinite.ind]
    A <- (8*nu.sub^2+12*nu.sub+3)/(96*(lambda.sub^(1/nu.sub))*nu.sub^2)
    B <- (1+6*nu.sub)/(144*(lambda.sub^(2/nu.sub))*nu.sub^3)
    sum1.sub <- nu.sub*lambda.sub^(1/nu.sub)-
      ((nu.sub-1)/(2*nu.sub)*log(lambda.sub)+(nu.sub-1)/2*log(2*pi)+0.5*log(nu.sub))+
      log(1+(nu.sub-1)*(A+B))
    sum1[infinite.ind] <- sum1.sub
  }
  if (!log.z){
    sum1 <- exp(sum1)}
  return(sum1)
}


logZ <- function(log_lambda, nu, summax = 100){
  # approximates normalizing constant for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(log_lambda=log_lambda, nu=nu)
  log_lambda <- df[,1]
  nu <- df[,2]
  return(logZ_c(log_lambda,nu,summax=summax))
}
