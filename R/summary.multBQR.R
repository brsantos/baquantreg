## Getting the results of each model

summary.multBQR <- function(object, ...){
  if (length(object$modelsDir) == 1){
    if (object$method == 'rcpp')
      lapply(object, summary.bqr, mult = TRUE, ...)
    else {
      lapply(object$modelsDir, function(a){
        lapply(a$modelsTau, function(b){
          dimfix <- dim(b$fixed.effects)
          BetaPosterior <- b$fixed.effects[1:dimfix[1], 1:dimfix[2]]

        })
      })
    }
  }
}