#' Two part model using Bayesian quantile regression model
#'
#' This function estimates a two part model using a Bayesian quantile
#' regression model to describe the continous part of the conditional
#' distribution. The response variable is assumed to follow a mixed
#' discrete-continuous distribution.
#'
#' @param formula a formula object, with the response on the left of a ~
#'  operator, and the terms, separated by + operators, on the right. If the
#'  assumed point mass distribution is at 1, then all values equal to 1 should
#'  be replaced by 0.
#' @param tau Quantile of interest.
#' @param data a data.frame from which to find the variables defined in the
#' formula
#' @param itNum Number of iterations.
#' @param thin Thinning parameter.
#' @param betaValue Initial values for the parameter beta for the continuous part.
#' @param sigmaValue Initial value for the scale parameter.
#' @param vSampleInit ?
#' @param gammaValue Initial value for the parameter gamma of the discrete
#'  part.
#' @param sigmaGamma Tuning parameter for the Metropolis-Hastings step.
#' @param link Integer defining the link function used for the probability
#'  model. Default is 1 for the logit link function. If 0, the probit function
#'  is used.
#' @param priorVar Value that multiplies a identity matrix in the elicition
#'  process of the prior variance of the regression parameters.
#' @param refresh Interval between printing a message during the iteration
#'  process. Default is set to 100.
#' @param quiet Logical. If FALSE it will print messages depending on the
#'  refresh parameter to show that the chain is updating. If TRUE it will not
#'  print messages during the iteration process.
#' @return A list with the chains of all parameters of interest.
#' @references Santos and Bolfarine (2015) - Bayesian quantile regression
#'  analysis for continuous data with a discrete component at zero.
#' \emph{Preprint}. \url{http://arxiv.org/abs/1511.05925}
#' @export
#' @useDynLib baquantreg
#' @examples
#' set.seed(1)
#' data('BrazilElectricity')
#' modelo2p <- twopartQR(prop_elec/100 ~ population + income_percap, tau=0.5,
#'  data=BrazilElectricity, itNum=2000, sigmaGamma=2, quiet = TRUE)

twopartQR <- function(formula, tau = 0.5, data, itNum, thin=1,
                      betaValue = NULL, sigmaValue=1,
                      vSampleInit = NULL, gammaValue = NULL,
                      sigmaGamma = 0.5, link=1, priorVar = 100,
                      refresh = 100, quiet = FALSE){

  y <- as.numeric(stats::model.extract(stats::model.frame(formula, data), 'response'))
  X <- stats::model.matrix(formula, data)

  if (is.null(betaValue)) betaValue <- rep(0, dim(X)[2])
  if (is.null(gammaValue)) gammaValue <- rep(0, dim(X)[2])
  if (is.null(vSampleInit)) vSampleInit <- rep(1, sum(y!=0))

  twoPartModel <- list()
  twoPartModel$chains <- lapply(tau, function(a){
    tpBayesQR(tau = a, y=y, X=X, itNum=itNum, thin=thin,
              betaValue=betaValue, sigmaValue=sigmaValue,
              vSampleInit = vSampleInit, gammaValue = gammaValue,
              sigmaGamma=sigmaGamma, link=link, priorVar = priorVar,
              refresh=refresh, quiet = quiet)
  })

  twoPartModel$tau <- tau
  twoPartModel$formula <- formula
  twoPartModel$data <- data
  twoPartModel$link <- link

  class(twoPartModel) <- "twopartQR"
  return(twoPartModel)
}
