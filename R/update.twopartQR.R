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
  if(class(model)!='twopartQR')
    stop("Class different of 'twopartQR'. Use a different update function.")

  newChains <- lapply(1:length(model$tau), function(a){
    modelFit <- twopartQR(model$formula, model$tau[a], model$data, itNum+1,
                        betaValue=tail(model$chains[[a]]$BetaSample,1),
                        sigmaValue=tail(model$chains[[a]]$SigmaSample, 1),
                        vSampleInit=tail(model$chains[[a]]$vSample, 1),
                        link = model$link, quiet = T)$chains

    list(BetaSample = modelFit[[1]]$BetaSample,
         SigmaSample = modelFit[[1]]$SigmaSample,
         GammaSample = modelFit[[1]]$GammaSample,
         vSample = modelFit[[1]]$vSample,
         acceptRate = modelFit[[1]]$acceptRate)
  })

  dimChains <- dim(model$chains[[1]]$BetaSample)[1]
  w1 <- dimChains/(itNum + dimChains)
  w2 <- 1 - w1

  newTwopart <- list()
  newTwopart$chains <- lapply(1:length(model$tau), function(aa){
      BetaSample = rbind(model$chains[[aa]]$BetaSample,
                         newChains[[aa]]$BetaSample[-1,])
      SigmaSample = c(model$chains[[aa]]$SigmaSample,
                      newChains[[aa]]$SigmaSample[-1])
      vSample = rbind(model$chains[[aa]]$vSample,
                      newChains[[aa]]$vSample[-1,])
      GammaSample = c(model$chains[[aa]]$GammaSample,
                      newChains[[aa]]$GammaSample[-1])
      acceptRate = w1*model$chains[[aa]]$acceptRate +
        w2*newChains[[aa]]$acceptRate

      list(BetaSample = BetaSample, SigmaSample = SigmaSample,
           vSample = vSample, GammaSample = GammaSample,
           acceptRate = acceptRate)
  })

  newTwopart$tau <- model$tau
  newTwopart$formula <- model$formula
  newTwopart$data <- model$data
  newTwopart$link <- model$link

  class(newTwopart) <- "twopartQR"
  return(newTwopart)
}