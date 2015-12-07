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
#' @seealso zitobitQR
#' @export
#' @S3method update zitobitQR
#' @examples
#' set.seed(1)
#'

update.zitobitQR <- function(model, numIt, thin, sigmaGamma){
  if(class(model)!='zitobitQR')
    stop("Class different of 'zitobitQR'. Use a different update function.")

  newChains <- zitobitQR(tau = model$tau, y=model$y, X=model$X, numIt=numIt,
                         thin=thin, betaValue=tail(model$betaSample,1),
                         sigmaValue=tail(model$sigmaSample, 1),
                         gammaValue=tail(model$gammaSample, 1),
                         sigmaGamma = sigmaGamma, link=model$link)

  w1 <- dim(model$betaSample)[1]/(numIt + dim(model$betaSample)[1])
  w2 <- numIt/(numIt + dim(model$betaSample)[1])

  newZiTobit <- list(tau = model$tau, y=model$y, X=model$X,
                     betaSample = rbind(model$betaSample,
                                        newChains$betaSample),
                     sigmaSample = rbind(model$sigmaSample,
                                         newChains$sigmaSample),
                     gammaSample = rbind(model$gammaSample,
                                         newChains$gammaSample),
                     acceptRate = w1*model$acceptRate + w2*newChains$acceptRate)
  class(newZiTobit) <- "zitobitQR"

  return(newZiTobit)
}