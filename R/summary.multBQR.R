## Getting the results of each model

summary.multBQR <- function(object, ...){
   lapply(object, summary.bqr, mult = TRUE, ...)
}