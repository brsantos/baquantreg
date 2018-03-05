## Getting the results of each model

summary.multBQR <- function(object, ...){
   lapply(object, baquantreg:::summary.bqr, mult = TRUE, ...)
}

