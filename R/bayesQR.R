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
#' @param priorVar Value that multiplies a identity matrix in the elicition
#' process of the prior variance of the regression parameters.
#' @param refresh Interval between printing a message during the iteration
#' process. Default is set to 100.
#' @param quiet If TRUE, the default, it does not print messages to check if
#'  the MCMC is actually updating. If FALSE, it will use the value of refresh
#'  to print messages to control the iteration process.
#' @param tobit If TRUE, it will input the censored value for all observations
#'  with y = 0, according to the model. If FALSE, the default, it will estimate
#'  the parameter without this inputation process.
#' @param recordLat If TRUE, it will keep the Markov chain samples for the
#'  latent variable. Default is FALSE.
#' @param blocksV Number of blocks used to sample in the posterior distribution
#'  of the latent variable. If 0, then blocking is not used and all latent
#'  observations are sampled from. Default value is 0.
#' @param stopOrdering If TRUE, it will stop ordering the weighted residuals
#'  in order to update the states of the latent variables, and will consider
#'  the ordering of some particular state of the chain; if FALSE, for every
#'  iteration of the MCMC procedure, it will keep reordering these residual
#'  terms. Default is FALSE.
#' @param numOrdered The number of iterations that will be used to order
#'  the weighted residuals needed for the update of the posterior
#'  distribution of the latent variables. Default is half the size of
#'  the MCMC chain.
#' @return A list with the chains of all parameters of interest.
#' @references Kozumi and Kobayashi (2011) - Gibbs sampling methods for
#'  Bayesian quantile regression. Journal of Statistical Computation and
#'  Simulation.
#' @export
#' @useDynLib baquantreg
#' @examples
#' set.seed(1)

bayesQR <- function(formula, tau = 0.5, data, itNum, thin=1,
                    betaValue = NULL, sigmaValue=1, vSampleInit = NULL,
                    priorVar = 100, refresh = 100, quiet = T,
                    tobit = FALSE, recordLat = FALSE, blocksV = 0,
                    stopOrdering = FALSE, numOrdered = itNum/2){

  y <- as.numeric(stats::model.extract(stats::model.frame(formula, data), 'response'))
  X <- stats::model.matrix(formula, data)

  if (is.null(betaValue)) betaValue <- rep(0, dim(X)[2])
  if (is.null(vSampleInit)) vSampleInit <- rep(1, length(y))

  output <- list()
  output$chains <- lapply(tau, function(a){
    BayesQR(tau = a, y = y, X = X, itNum = itNum, thin = thin,
            betaValue = betaValue, sigmaValue = sigmaValue,
            vSampleInit = vSampleInit, priorVar = priorVar,
            refresh = refresh, quiet = quiet, tobit = tobit,
            recordLat = recordLat, blocksV = blocksV)
  })

  output$tau <- tau
  output$formula <- formula
  output$data <- data

  class(output) <- "bqr"

  return(output)
}
