#' Model update for a Bayesian quantile regression model
#'
#' This function updates the chains for a model fitted by the function bayesQR.
#'
#' @param model Fitted model by the function bayesQR
#' @param numIt Size of the chain that will be added to the current available
#' chain.
#' @param thin Thinning parameter for the chains.
#' @return A list with the updated chains.
#' @seealso \code{\link{twopartQR}}
#' @export
#' @examples
#' set.seed(1)
#'

update.bqr <- function(model, itNum, ...){
  if(class(model)!='bqr')
    stop("Class different from 'bqr'. Use a different update function.")

  newChains <- lapply(1:length(model$tau), function(a){
    modelFit <- bayesQR(model$formula, model$tau[a], model$data, itNum+1,
            betaValue=tail(model$chains[[a]]$BetaSample,1),
            sigmaValue=tail(model$chains[[a]]$SigmaSample, 1),
            vSampleInit=tail(model$chains[[a]]$vSample, 1))$chains

    list(BetaSample = modelFit[[1]]$BetaSample,
         SigmaSample = modelFit[[1]]$SigmaSample,
         vSample = modelFit[[1]]$vSample)
  })

  newBayesQR <- list()
  newBayesQR$chains <- lapply(1:length(model$tau), function(aa){
    BetaSample = rbind(model$chains[[aa]]$BetaSample,
                       newChains[[aa]]$BetaSample[-1,])
    SigmaSample = c(model$chains[[aa]]$SigmaSample,
                        newChains[[aa]]$SigmaSample[-1])
    vSample = rbind(model$chains[[aa]]$vSample,
                    newChains[[aa]]$vSample[-1,])

   list(BetaSample = BetaSample, SigmaSample = SigmaSample,
        vSample = vSample)
  })

  newBayesQR$formula <- model$formula
  newBayesQR$tau = model$tau
  newBayesQR$data = model$data

  class(newBayesQR) <- "bqr"
  return(newBayesQR)
}