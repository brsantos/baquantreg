#' Probability of being an outlier given the Bayesian quantile regression fits.
#'
#' Returns the probability that each observation is an outlier given the
#'  Bayesian quantile regression fits.
#'
#' @param object This is an object of the class "bqr", produced by a call
#'  to the bayesQR function.
#' @param burnin Initial part of the chain, which is to be discarded. Default
#'  value is 50.
#' @param plotProb If TRUE, the function prints the plot with all probabilities.
#'  Default is set to TRUE.
#' @param scales.free If FALSE, the default, then all plots will use the same
#'  y scale. If TRUE, for each tau the plot will use the best possible scale
#'  in order to visualize the probability information for all observations.
#' @return Prints a plot of the posterior probability of being an outlier for
#'  all observations and returns a data.frame with all values of the
#'  probabilities.
#' @export
#' @useDynLib baquantreg
#' @import ggplot2
#' @importFrom FNN KL.divergence

latentKL <- function(object, burnin = 50, plotKL = TRUE,
                        scales.free = FALSE){

  if (class(object) != "bqr")
    stop("This function is not suited for your model.")

  taus <- object$tau
  nobs <- dim(object$chains[[1]]$vSample)[2]

  klValues <- t(sapply(1:nobs, function(b){
    as.numeric(sapply(object$chains, function(a){
      sapply(1:nobs, function(c){
        if (b != c){
          max(FNN::KL.divergence(a$vSample[-c(1:burnin), b],
                                 a$vSample[-c(1:burnin), c], k = 5))
        }
        else 0
      })
    }))
  }))

  plotData <- data.frame(nobs = rep(1:nobs, times=dim(klValues)[2]),
                         values = as.numeric(klValues),
                         taus = rep(taus, each=dim(klValues)[1]))

  maxProb <- which.max(aggregate(values ~ nobs, data=plotData, mean)$values)
  print(paste("The observation with greater Kullback-Leibler divergence
              from the others is:",
              maxProb))

  g <- ggplot(subset(plotData), aes(y=values, x=nobs)) + theme_bw()
  if (length(taus) > 1){
    if (scales.free) g <- g + facet_wrap(~taus, scales='free')
    else g <- g + facet_wrap(~taus)
  }

  g <- g + geom_point() +
    ylab("Mean Kullback-Leibler divergence for each point") + xlab("# Observation")

  if (plotKL) print(g)

  return(plotData)
}
