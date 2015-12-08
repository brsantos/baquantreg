#' Model update
#'
#' This function updates the chains for a model fitted by the function
#' zitobitQR.
#'
#' @param model Fitted model by the function twopartQR
#' @param numIt Size of the chain that will be added to the current available
#' chain.
#' @param thin Thinning parameter for the chains.
#' @param sigmaGamma Tuning parameter for the Metropolis-Hasting part of the
#' algorithm.
#' @return A list with the updated chains.
#' @seealso \code{\link{twopartQR}}
#' @export
#' @examples
#' set.seed(1)
#'

update.twopartQR <- function(model, itNum, thin, sigmaGamma, ...){
  if(class(model)!='zitobitQR')
    stop("Class different of 'twopartQR'. Use a different update function.")

  newChains <- twopartQR(model$formula, model$tau, model$data, itNum,
                         thin, betaValue=tail(model$betaSample,1),
                         sigmaValue=tail(model$sigmaSample, 1),
                         gammaValue=tail(model$gammaSample, 1),
                         sigmaGamma = sigmaGamma, link=model$link, ...)

  w1 <- dim(model$betaSample)[1]/(itNum + dim(model$betaSample)[1])
  w2 <- itNum/(itNum + dim(model$betaSample)[1])

  newTwopart <- list(tau = model$tau, y=model$y, X=model$X,
                     betaSample = rbind(model$betaSample,
                                        newChains$betaSample),
                     sigmaSample = rbind(model$sigmaSample,
                                         newChains$sigmaSample),
                     gammaSample = rbind(model$gammaSample,
                                         newChains$gammaSample),
                     acceptRate = w1*model$acceptRate +
                       w2*newChains$acceptRate)
  class(newTwopart) <- "twopartQR"

  return(newTwopart)
}