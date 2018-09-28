#' Combine R Objects by Columns
#' 
#' Take a sequence of vector, matrix or data-frame arguments and combine them by columns. 
#' \code{CBIND} is used within the package over \code{cbind} to recycle the shorter 
#' arguments so that their number of rows would match. 
#' 
#' @param ... (generalized) vectors or matrices. These can be given as named arugments. 
#' @param deparse.level integer; deparse.level = 0 constructs no labels, 
#' deparse.level = 1 (the default) or > 1 constructs labels from the arguments names. 
#' 
CBIND <- function(..., deparse.level = 1) {
  dots <- list(...)
  len <- sapply(dots, length)
  dots <- lapply(seq_along(dots),
                 function(i, x, len) rep(x[[i]], length.out = len),
                 x = dots, len = max(len))
  do.call(cbind, c(dots, deparse.level = deparse.level))
}

#' Test for a whole number
#' 
#' Test for integer/whole number vector 
#' 
#' @param x numeric vector to be tested  
#' @param tol numeric; precision level
#' 
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}

