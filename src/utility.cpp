#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the Normalizing Constant in log scale for COM-Poisson distribution
//'  
//' The calculation of the function \code{logZ} will be performed here. This function is 
//' used to approximate the normalizing constant for COM-Poisson distributions via 
//' truncation. The standard COM-Poisson parametrization is being used here. 
//' 
//' It is assumed that vectors log_lambda & nu are of equal length. 
//' 
//' This function was originally purposed in the \code{cmpreg} package of Ribeiro Jr, 
//' Zeviani & Demétrio (2019).
//' 
//' @param log_lambda rate parameter in log scale.
//' @param nu dispersion parameter, straightly positive.
//' @param summax maximum number of terms to be considered in the truncated sum.
//' @useDynLib mpcmp, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @references 
//' Ribeiro Jr, E. E., Zeviani, W. M., Demétrio, C. G. B. (2019) \code{cmpreg}: 
//' Reparametrized COM-Poisson Regression Models. R package version 0.0.1.
//' @export
//' 
// [[Rcpp::export]]
NumericVector logZ_c(NumericVector log_lambda, NumericVector nu, int summax) {
  // Control loop
  // int maxiter = 1e4;
  double log_epsilon = std::log(1e-10);
  // Output vector
  int n = log_lambda.size();
  NumericVector out(n);
  // Compute logz
  for (int i = 0; i < n; ++i) {
    double logz  = 0;
    double logz_ = 0;
    for (int j = 1; j < summax; ++j) {
      logz_ += log_lambda[i] - nu[i] * log(j);
      logz = R::logspace_add(logz, logz_);
      if (logz_ - logz < log_epsilon) break;
    }
    out[i] = logz;
  }
  return out;
}
