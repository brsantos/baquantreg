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
#' @param vSampleInit Initial value for the latent variables.
#' @param priorVar Value that multiplies a identity matrix in the elicitation
#'  process of the prior variance of the regression parameters.
#' @param hyperSigma Vector of size containing the hyperparameters of the
#'  inverse gamma distribution for the sigma parameter of the asymmetric
#'  Laplace distribution. Default is c(0.1, 0.1), which gives a noninformative
#'  prior for sigma.
#' @param refresh Interval between printing a message during the iteration
#' process. Default is set to 100.
#' @param sigmaSampling If TRUE, the default, the MCMC procedure will draw values
#'  from the posterior distribution of sigma. Otherwise, it will fix the
#'  value to 1 for all values of the chain.
#' @param quiet If TRUE, the default, it does not print messages to check if
#'  the MCMC is actually updating. If FALSE, it will use the value of refresh
#'  to print messages to control the iteration process.
#' @param tobit If TRUE, it will input the censored value for all observations
#'  with y = 0, according to the model. If FALSE, the default, it will estimate
#'  the parameter without this inputation process.
#' @param recordLat If TRUE, it will keep the Markov chain samples for the
#'  latent variable. Default is FALSE.
#' @return A list with the chains of all parameters of interest.
#' @references Kozumi and Kobayashi (2011) - Gibbs sampling methods for
#'  Bayesian quantile regression. Journal of Statistical Computation and
#'  Simulation.
#' @export
#' @useDynLib baquantreg
#' @examples
#' set.seed(1)

bayesQR <- function(formula, tau = 0.5, data, itNum, thin=1,
                    betaValue = NULL, sigmaValue = 1, vSampleInit = NULL,
                    priorVar = 100, hyperSigma = c(0.1, 0.1),
                    refresh = 100, sigmaSampling = TRUE,
                    quiet = T,
                    tobit = FALSE, recordLat = FALSE){

  y <- as.numeric(
    stats::model.extract(stats::model.frame(formula, data), 'response')
    )
  X <- stats::model.matrix(formula, data)

  if (is.null(betaValue)) betaValue <- rep(0, dim(X)[2])
  if (is.null(vSampleInit)) vSampleInit <- rep(1, length(y))

  output <- list()
  output$chains <- lapply(tau, function(a){
    BayesQR(tau = a, y = y, X = X, itNum = itNum, thin = thin,
            betaValue = betaValue, sigmaValue = sigmaValue,
            vSampleInit = vSampleInit, priorVar = priorVar,
            hyperSigma = hyperSigma,
            refresh = refresh, sigmaSampling = sigmaSampling,
            quiet = quiet, tobit = tobit,
            recordLat = recordLat)
  })

  output$tau <- tau
  output$formula <- formula
  output$data <- data

  class(output) <- "bqr"

  return(output)
}
