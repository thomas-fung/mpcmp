#' Calculate the Normalizing Constant for COM-Poisson distribution
#' 
#' A function to approximate the normalizing constant for COM-Poisson distributions
#' via truncation. The standard COM-Poisson parametrization is being used here. 
#' Notice that the sum is hard coded to tuncate at 100 so the approximation will be quite
#' bad if the COM-Poisson has a large rate or mean. 
#' 
#' @param lambda rate parameter, straightly positive
#' @param nu diepsersoin parameter, straightly positive
Z <- function(lambda, nu){
  # approximates normalizing constant for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(lambda=lambda, nu=nu)
  lambda <- df[,1]
  nu <- df[,2]
  summax <- 100
  termlim <- 1e-6
  # zero vector length of same length as lambda
  sum1 <- 0
  for(y in 1:summax) {
    term <- exp(log(lambda^(y-1)) - nu*lgamma(y))
    if (y > 3){
      if (max(term/sum1,na.rm = TRUE) < termlim) {
        break
      }
    }
    sum1 <- sum1 + term
  }
  return(sum1)
}


