#include <Rcpp.h>
using namespace Rcpp;

//' @title Compute the normalising constant
//' @description Computes the normalising constant
//' @details More description here
//' @param log_lambda: 
//' @export
// [[Rcpp::export]]
NumericVector compute_log_constant(NumericVector log_lambda, NumericVector nu, int maxiter) {
  // Control loop
  // int maxiter = 1e4;
  double log_epsilon = std::log(1e-8);
  // Output vector
  int n = log_lambda.size();
  NumericVector out(n);
  // Compute logz
  for (int i = 0; i < n; ++i) {
    double logz  = 0;
    double logz_ = 0;
    for (int j = 1; j < maxiter; ++j) {
      logz_ += log_lambda[i] - nu[i] * log(j);
      logz = R::logspace_add(logz, logz_);
      if (logz_ - logz < log_epsilon) break;
    }
    out[i] = logz;
  }
  return out;
}