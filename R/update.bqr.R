#' Model update for a Bayesian quantile regression model
#'
#' This function updates the chains for a model fitted by the function bayesQR.
#'
#' @param object Fitted model by the function bayesQR
#' @param itNum Size of the chain that will be added to the current available
#' chain.
#' @param ... other summary parameters (currently not used)
#' @return A list with the updated chains.
#' @seealso \code{\link{twopartQR}}
#' @export
#' @examples
#' set.seed(1)

update.bqr <- function(object, itNum, ...){
  if(class(object)!='bqr')
    stop("Class different from 'bqr'. Use a different update function.")

  newChains <- lapply(1:length(object$tau), function(a){
    modelFit <- bayesQR(object$formula, object$tau[a], object$data, itNum+1,
            betaValue=tail(object$chains[[a]]$BetaSample,1),
            sigmaValue=tail(object$chains[[a]]$SigmaSample, 1),
            vSampleInit=tail(object$chains[[a]]$vSample, 1))$chains

    list(BetaSample = modelFit[[1]]$BetaSample,
         SigmaSample = modelFit[[1]]$SigmaSample,
         vSample = modelFit[[1]]$vSample)
  })

  newBayesQR <- list()
  newBayesQR$chains <- lapply(1:length(object$tau), function(aa){
    BetaSample = rbind(object$chains[[aa]]$BetaSample,
                       newChains[[aa]]$BetaSample[-1,])
    SigmaSample = c(object$chains[[aa]]$SigmaSample,
                        newChains[[aa]]$SigmaSample[-1])
    vSample = rbind(object$chains[[aa]]$vSample,
                    newChains[[aa]]$vSample[-1,])

   list(BetaSample = BetaSample, SigmaSample = SigmaSample,
        vSample = vSample)
  })

  newBayesQR$formula <- object$formula
  newBayesQR$tau = object$tau
  newBayesQR$data = object$data

  class(newBayesQR) <- "bqr"
  return(newBayesQR)
}