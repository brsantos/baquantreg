#' Two part model using Bayesian quantile regression model
#'
#' This function estimates a two part model using a Bayesian quantile
#' regression model to describe the continous part of the conditional
#' distribution. The response variable is assumed to follow a mixed
#' discrete-continuous distribution.
#'
#' @param formula a formula object, with the response on the left of a ~
#'  operator, and the terms, separated by + operators, on the right.
#' @param tau Quantile of interest.
#' @param data a data.frame from which to find the variables defined in the
#' formula
#' @param itNum Number of iterations.
#' @param thin Thinning parameter.
#' @param betaValue Initial values for the parameter beta for the continuous part.
#' @param sigmaValue Initial value for the scale parameter.
#' @param gammaValue Initial value for the parameter gamma of the discrete
#' part.
#' @param sigmaGamma Tuning parameter for the Metropolis-Hastings step.
#' @param link Integer defining the link function used for the probability
#' model. Default is 1 for the logit link function.
#' @param priorVar Value that multiplies a identity matrix in the elicition
#' process of the prior variance of the regression parameters.
#' @param refresh Interval between printing a message during the iteration
#' process. Default is set to 100.
#' @return A list with the chains of all parameters of interest.
#' @references Santos and Bolfarine (2015) - Bayesian quantile regression
#' analysis for continuous data with a discrete component at zero.
#' \emph{Preprint}. \url{http://arxiv.org/abs/1511.05925}
#' @export
#' @useDynLib baquantreg
#' @examples
#' set.seed(1)

twopartQR <- function(formula, tau = 0.5, data, itNum, thin=1,
                      betaValue = NULL, sigmaValue=1, gammaValue = NULL,
                      sigmaGamma = 0.5, link=1, priorVar = 100,
                      refresh = 100){

  y <- as.numeric(model.extract(model.frame(formula, data), 'response'))
  X <- model.matrix(formula, data)

  if (is.null(betaValue)) betaValue <- rep(0, dim(X)[2])
  if (is.null(gammaValue)) gammaValue <- rep(0, dim(X)[2])

  twoPartModel <- tpBayesQR(tau = tau, y=y, X=X, itNum=itNum, thin=thin,
                          betaValue=betaValue, sigmaValue=sigmaValue,
                          gammaValue = gammaValue, sigmaGamma=sigmaGamma,
                          link=link, priorVar = priorVar, refresh=refresh)

  twoPartModel$tau <- tau
  twoPartModel$formula <- formula
  twoPartModel$data <- data
  twoPartModel$link <- link

  class(twoPartModel) <- "twopartQR"

  return(twoPartModel)
}
