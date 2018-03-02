#' Summary methods for Bayesian quantile regression models.
#'
#' Returns a summary data.frame for a Bayesian quantile regression fit for more
#'  than one quantile.
#'
#' @param object This is an object of class "bqr", produced by a call to the
#'  bayesQR function.
#' @param burnin Initial part of the chain, which is to be discarded. Default
#'  value is 1000.
#' @param ci Credible interval coefficient. Default value is 0.95.
#' @param mult Defining whether this function will be used for multiple-output
#'  Bayesian quantile regression or not. Default is FALSE.
#' @param ... other summary parameters (not used)
#' @return A data frame with summary information about the quantile regression
#'  parameters.
#' @export
#' @useDynLib baquantreg
#' @examples
#' set.seed(1)


summary.bqr <- function (object, burnin = 1000, ci = 0.95, mult = FALSE, ...)
{
  if (class(object) != "bqr")
    stop("Use the correct summary method for your model")

  numIt <- dim(object$chains[[1]]$BetaSample)[1]

  X <- model.matrix(object$formula, object$data)

  output <- list()
  output$BetaPosterior <- lapply(object$chains, function(a){

    if (!mult) vnames <- colnames(X)
    else vnames <- c(colnames(X), 'directionX')

    coef <- apply(a$BetaSample[(burnin+1):numIt, ], 2, mean)
    quantilesL <- apply(a$BetaSample[(burnin+1):numIt, ], 2, quantile, (1-ci)/2)
    quantilesU <- apply(a$BetaSample[(burnin+1):numIt, ], 2, quantile, 1-(1-ci)/2)

    data.frame(variable = vnames,
               coef = coef,
               lower = quantilesL,
               upper = quantilesU)
  })

  names(output$BetaPosterior) <- paste("Tau = ", object$tau)

  output$SigmaPosterior <- data.frame(t(sapply(object$chains, function(a){
    meanSigma <- mean(a$SigmaSample[(burnin+1):numIt])
    quantilesL <- quantile(a$SigmaSample[(burnin+1):numIt], (1-ci)/2)
    quantilesU <- quantile(a$SigmaSample[(burnin+1):numIt], 1-(1-ci)/2)

    c(meanSigma, quantilesL, quantilesU)
  })))

  colnames(output$SigmaPosterior) <- c("Mean", "Lower", "Upper")

  tau <- object$tau
  output$SigmaPosterior <- data.frame(cbind(tau, output$SigmaPosterior))

  output$taus <- object$tau

  class(output) <- "summary.bqr"
  output
}