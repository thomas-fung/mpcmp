#' @title Compute the normalising constant
#' @description Computes the normalising constant
#' @export
compute_log_constant <- function(log_lambda, nu, maxiter = 100) {
    .Call('_mpcmp_compute_log_constant', PACKAGE = 'mpcmp', log_lambda, nu, maxiter)
}

