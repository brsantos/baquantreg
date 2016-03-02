#' Bayesian quantile regression model with a discrete component at zero
#'
#' This function estimates a bayesian quantile regression model with a discrete component
#' at zero, where all zero observations are assumed to distributed according to a mixed
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
#' @param gamaValue Initial value for the parameter gamma of the discrete part.
#' @param sigmaBetaZero Tuning parameter for the Metropolis-Hastings step.
#' @param link Integer defining the link function used for the probability
#' model. Default is 1.
#' for the logit link function.
#' @param priorVar Value that multiplies a identity matrix in the elicition
#' process of the prior variance of the regression parameters.
#' @param refresh Interval between printing a message during the iteration
#' process. Default is set to 100.
#' @param quiet Logical. If FALSE (default) it will print messages depending on
#'  the refresh parameter to show that the Markov chain is updating. If TRUE it
#'  will not print messages during the iteration process.
#' @return A list with the chains of all parameters of interest.
#' @references Santos and Bolfarine (2015) - Bayesian quantile regression
#'  analysis for continuous data
#' with a discrete component at zero. \emph{Preprint}.
#'  \url{http://arxiv.org/abs/1511.05925}
#' @export
#' @useDynLib baquantreg
#' @import RcppArmadillo
#' @import RcppGSL
#' @importFrom Rcpp evalCpp
#' @examples
#' set.seed(1)
#' data("BrazilDurableGoods")
#' # Change the number of iterations for better results.
#' model <- zitobitQR(expenditure ~ age + education, tau=0.5,
#'                   data=BrazilDurableGoods, itNum=100,
#'                   sigmaGamma=0.10, refresh=20)

zitobitQR <- function(formula, tau = 0.5, data, itNum, thin=1,
                      betaValue = NULL, sigmaValue=1,
                      gammaValue = NULL,  sigmaGamma = 0.5,
                      link=1, priorVar = 100, refresh = 100,
                      quiet = FALSE){

  y <- as.numeric(model.extract(model.frame(formula, data), 'response'))
  X <- model.matrix(formula, data)

  if (is.null(betaValue)) betaValue <- rep(0, dim(X)[2])
  if (is.null(gammaValue)) gammaValue <- rep(0, dim(X)[2])

  ziTobit <- list()

  ziTobit$chains <- ziTobitBayesQR(tau = tau, y=y, X=X, itNum=itNum, thin=thin,
                          betaValue=betaValue, sigmaValue=sigmaValue,
                          gammaValue=gammaValue, sigmaGamma=sigmaGamma,
                          link=link, priorVar=priorVar, refresh=refresh,
                          quiet=quiet)

  ziTobit$tau <- tau
  ziTobit$formula <- formula
  ziTobit$data <- data
  ziTobit$link <- link

  class(ziTobit) <- "zitobitQR"

  return(ziTobit)
}
