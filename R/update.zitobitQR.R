#' Model update
#'
#' This function updates the chains for a model fitted by the function
#' zitobitQR.
#'
#' @param object Fitted model by the function zitobitQR
#' @param itNum Size of the chain that will be added to the current available
#'   chain.
#' @param thin Thinning parameter for the chains.
#' @param sigmaGamma Tuning parameter for the Metropolis-Hasting part of the
#'   algorithm.
#' @param ... other summary parameters (currently not used)
#' @return A list with the updated chains.
#' @seealso \code{\link{zitobitQR}}
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' data("BrazilDurableGoods")
#' # Change the number of iterations for better results.
#' model <- zitobitQR(expenditure ~ age + education, tau = 0.5,
#'                   data = BrazilDurableGoods, itNum = 100,
#'                   sigmaGamma = 0.10, refresh = 20)
#' model <- update(model, 100, thin=1, sigmaGamma = 0.10, refresh=10)
#' }
update.zitobitQR <- function(object, itNum, thin=1, sigmaGamma, ...){
  if(class(model)!='zitobitQR')
    stop("Class different of 'zitobitQR'. Use a different update function.")

  newChains <- zitobitQR(object$formula, object$tau, object$data, itNum=itNum,
                         thin=thin, betaValue=tail(object$betaSample,1),
                         sigmaValue=tail(object$sigmaSample, 1),
                         gammaValue=tail(object$gammaSample, 1),
                         sigmaGamma = sigmaGamma, link=object$link, ...)

  w1 <- dim(object$betaSample)[1]/(itNum + dim(object$betaSample)[1])
  w2 <- itNum/(itNum + dim(object$betaSample)[1])

  newZiTobit <- list(formula = object$formula, tau = object$tau,
                     betaSample = rbind(object$betaSample,
                                        newChains$betaSample),
                     sigmaSample = rbind(object$sigmaSample,
                                         newChains$sigmaSample),
                     gammaSample = rbind(object$gammaSample,
                                         newChains$gammaSample),
                     acceptRate = w1*object$acceptRate +
                       w2*newChains$acceptRate)
  class(newZiTobit) <- "zitobitQR"

  return(newZiTobit)
}