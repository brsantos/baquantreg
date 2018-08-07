## Getting the results of each model

summary.multBQR <- function(object, ...){
  if (object$method == 'rcpp')
    lapply(object, summary.bqr, mult = TRUE, ...)
  else {

  }
}