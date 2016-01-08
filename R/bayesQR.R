#' Bayesian quantile regression model
#'
#' This function estimates a Bayesian quantile regression model
#' The response variable is assumed to follow a asymmetric Laplace
#' distribution.
#'
#' @param formula a formula object, with the response on the left of a ~
#'  operator, and the terms, separated by + operators, on the right.
#' @param tau Quantile of interest.
#' @param data A data.frame from which to find the variables defined in the
#' formula
#' @param itNum Number of iterations.
#' @param thin Thinning parameter. Default value is 1.
#' @param betaValue Initial values for the parameter beta for the continuous
#' part.
#' @param sigmaValue Initial value for the scale parameter.
#' @param gammaValue Initial value for the parameter gamma of the discrete
#' part.
#' @param priorVar Value that multiplies a identity matrix in the elicition
#' process of the prior variance of the regression parameters.
#' @param refresh Interval between printing a message during the iteration
#' process. Default is set to 100.
#' @return A list with the chains of all parameters of interest.
#' @references Kozumi and Kobayashi (2015) - Gibbs sampling methods for
#' Bayesian quantile regression.
#' @export
#' @useDynLib baquantreg
#' @examples
#' set.seed(1)

bayesQR <- function(formula, tau = 0.5, data, itNum, thin=1,
                      betaValue = NULL, sigmaValue=1, priorVar = 100,
                      refresh = 100){

  y <- as.numeric(model.extract(model.frame(formula, data), 'response'))
  X <- model.matrix(formula, data)

  if (is.null(betaValue)) betaValue <- rep(0, dim(X)[2])

  output <- list()

  output$chains <- lapply(tau, function(a){
    BayesQR(tau = a, y=y, X=X, itNum=itNum, thin=thin,
            betaValue=betaValue, sigmaValue=sigmaValue,
            priorVar = priorVar,  refresh=refresh)
  })

  output$tau <- tau
  output$formula <- formula
  output$data <- data

  class(output) <- "bqr"

  return(output)
}
