#' Model update
#'
#' This function updates the chains for a model fitted by the function
#' zitobitQR.
#'
#' @param object Fitted model by the function twopartQR
#' @param itNum Size of the chain that will be added to the current available
#' chain.
#' @param thin Thinning parameter for the chains.
#' @param sigmaGamma Tuning parameter for the Metropolis-Hasting part of the
#' algorithm.
#' @param ... other summary parameters (currently not used)
#' @return A list with the updated chains.
#' @seealso \code{\link{twopartQR}}
#' @export
#' @examples
#' set.seed(1)
#'

update.twopartQR <- function(object, itNum, thin, sigmaGamma, ...){
  if(class(object)!='twopartQR')
    stop("Class different of 'twopartQR'. Use a different update function.")

  newChains <- lapply(1:length(object$tau), function(a){
    modelFit <- twopartQR(object$formula, object$tau[a], object$data, itNum+1,
                        betaValue=tail(object$chains[[a]]$BetaSample,1),
                        sigmaValue=tail(object$chains[[a]]$SigmaSample, 1),
                        vSampleInit=tail(object$chains[[a]]$vSample, 1),
                        link = object$link, quiet = T)$chains

    list(BetaSample = modelFit[[1]]$BetaSample,
         SigmaSample = modelFit[[1]]$SigmaSample,
         GammaSample = modelFit[[1]]$GammaSample,
         vSample = modelFit[[1]]$vSample,
         acceptRate = modelFit[[1]]$acceptRate)
  })

  dimChains <- dim(object$chains[[1]]$BetaSample)[1]
  w1 <- dimChains/(itNum + dimChains)
  w2 <- 1 - w1

  newTwopart <- list()
  newTwopart$chains <- lapply(1:length(object$tau), function(aa){
      BetaSample = rbind(object$chains[[aa]]$BetaSample,
                         newChains[[aa]]$BetaSample[-1,])
      SigmaSample = c(object$chains[[aa]]$SigmaSample,
                      newChains[[aa]]$SigmaSample[-1])
      vSample = rbind(object$chains[[aa]]$vSample,
                      newChains[[aa]]$vSample[-1,])
      GammaSample = c(object$chains[[aa]]$GammaSample,
                      newChains[[aa]]$GammaSample[-1])
      acceptRate = w1*object$chains[[aa]]$acceptRate +
        w2*newChains[[aa]]$acceptRate

      list(BetaSample = BetaSample, SigmaSample = SigmaSample,
           vSample = vSample, GammaSample = GammaSample,
           acceptRate = acceptRate)
  })

  newTwopart$tau <- object$tau
  newTwopart$formula <- object$formula
  newTwopart$data <- object$data
  newTwopart$link <- object$link

  class(newTwopart) <- "twopartQR"
  return(newTwopart)
}