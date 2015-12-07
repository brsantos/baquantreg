#' Bayesian quantile regression model with a discrete component at zero
#'
#' This function estimates a bayesian quantile regression model with a discrete component
#' at zero, where all zero observations are assumed to distributed according to a mixed
#' discrete-continuous distribution.
#'
#' @param tau Quantile of interest.
#' @param y Vector of response variable.
#' @param X Design matrix.
#' @param itNum Number of iterations.
#' @param thin Thinning parameter.
#' @param betaValue Initial values for the parameter beta for the continuous part.
#' @param sigmaValue Initial value for the scale parameter.
#' @param gamaValue Initial value for the parameter gamma of the discrete part.
#' @param sigmaBetaZero Tuning parameter for the Metropolis-Hastings step.
#' @param link Integer defining the link function used for the probability
#' model. Default is 1.
#' for the logit link function.
#' @return A list with the chains of all parameters of interest.
#' @references Santos and Bolfarine (2015) - Bayesian quantile regression analysis for continuous data
#' with a discrete component at zero. \emph{Preprint}. \url{http://arxiv.org/abs/1511.05925}
#' @export
#' @useDynLib baquantreg
#' @examples
#' set.seed(1)

zitobitQR <- function(tau = 0.5, y, X, itNum, thin=1,
                      betaValue=rep(0, dim(X)[2]),
                      sigmaValue=1, gammaValue=rep(0, dim(X)[2]),
                      sigmaGamma = 0.5, link=1){

  ziTobit <- ziTobitBayesQR(tau = tau, y=y, X=X, itNum=itNum, thin=thin,
                          betaValue=betaValue, sigmaValue=sigmaValue,
                          gammaValue = gammaValue,
                          sigmaGamma=sigmaGamma, link=link)

  ziTobit$tau <- tau
  ziTobit$y <- y
  ziTobit$X <- X
  ziTobit$link <- link

  class(ziTobit) <- "zitobitQR"

  return(ziTobit)
}


