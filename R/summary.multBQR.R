## Getting the results of each model

summary.multBQR <- function(object, ...){
  if (object$method == 'rcpp')
    lapply(object$modelsDir, summary.bqr, mult = TRUE, ...)
  else {
    lapply(object$modelsDir, function(a){
      lapply(a$modelsTau, function(b){
        colnamesEff <- colnames(b$fixed.effects)
        BetaPosterior <- b$fixed.effects[, colnamesEff]
        SigmaPosterior <- b$variance[, colnamesEff]
        list(BetaPosterior  = BetaPosterior, SigmaPosterior = SigmaPosterior)
      })
    })
  }
}