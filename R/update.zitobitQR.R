#' Model update
#'
#' This function updates the chains for a model fitted by the function
#' zitobitQR.
#'
#' @param model Fitted model by the function zitobitQR
#' @param numIt Size of the chain that will be added to the current available
#' chain.
#' @param thin Thinning parameter for the chains.
#' @param sigmaGamma Tuning parameter for the Metropolis-Hasting part of the
#' algorithm.
#' @return A list with the updated chains.
#' @seealso \code{\link{zitobitQR}}
#' @export
#' @examples
#' set.seed(1)
#' data("BrazilDurableGoods")
#' # Change the number of iterations for better results.
#' model <- zitobitQR(expenditure ~ age + education, tau=0.5,
#'                   data=BrazilDurableGoods, itNum=100,
#'                   sigmaGamma=0.10, refresh=20)
#' model <- update(model, 100, thin=1, sigmaGamma = 0.10, refresh=10)

update.zitobitQR <- function(model, itNum, thin=1, sigmaGamma, ...){
  if(class(model)!='zitobitQR')
    stop("Class different of 'zitobitQR'. Use a different update function.")

  newChains <- zitobitQR(model$formula, model$tau, model$data, itNum=itNum,
                         thin=thin, betaValue=tail(model$betaSample,1),
                         sigmaValue=tail(model$sigmaSample, 1),
                         gammaValue=tail(model$gammaSample, 1),
                         sigmaGamma = sigmaGamma, link=model$link, ...)

  w1 <- dim(model$betaSample)[1]/(itNum + dim(model$betaSample)[1])
  w2 <- itNum/(itNum + dim(model$betaSample)[1])

  newZiTobit <- list(formula = model$formula, tau = model$tau,
                     betaSample = rbind(model$betaSample,
                                        newChains$betaSample),
                     sigmaSample = rbind(model$sigmaSample,
                                         newChains$sigmaSample),
                     gammaSample = rbind(model$gammaSample,
                                         newChains$gammaSample),
                     acceptRate = w1*model$acceptRate +
                       w2*newChains$acceptRate)
  class(newZiTobit) <- "zitobitQR"

  return(newZiTobit)
}