#' Probability of being an outlier given the Bayesian quantile regression fits.
#'
#' Returns the probability that each observation is an outlier given the
#'  Bayesian quantile regression fits.
#'
#' @param object This is an object of the class "bqr", produced by a call
#'  to the bayesQR function.
#' @param burnin Initial part of the chain, which is to be discarded.
#'  Default value is 50.
#' @param plotProb If TRUE, the function prints the plot with all probabilities.
#'  Default is set to TRUE.
#' @param scales.free If FALSE, the default, then all plots will use the same
#'  y scale. If TRUE, for each tau the plot will use the best possible scale
#'  in order to visualize the probability information for all observations.
#' @param all.obs if TRUE, calculates KL divergence for all observations
#' @param obs if \code{all.obs} is FALSE, specifies the observation to
#'   calculate the KL divergence.
#' @return Prints a plot of the posterior probability of being an outlier for
#'  all observations and returns a data.frame with all values of the
#'  probabilities.
#' @export
#' @useDynLib baquantreg
#' @import ggplot2

probOutlier <- function(object, burnin = 50, plotProb = TRUE,
                        scales.free = FALSE, all.obs = TRUE,
                        obs){

  if (class(object) != "bqr")
    stop("This function is not suited for your model.")

  taus <- object$tau
  nobs <- dim(object$chains[[1]]$vSample)[2]

  if (all.obs) seqObs <- 1:nobs
  else seqObs <- obs

  prob <- sapply(object$chains, function(a){
    sapply(seqObs, function(b){
      maxValues <- apply(a$vSample[-c(1:burnin),-b], 2, max)
      vSample <- a$vSample[-c(1:burnin), b]
      mean(sapply(maxValues, function(c){
        sum(vSample > c)/length(vSample)
      }))
    })
  })

  plotData <- data.frame(nobs = rep(seqObs, times=length(taus)),
                           values = as.numeric(prob),
                           taus = rep(taus, each=length(seqObs)))

  if (all.obs){
    maxProb <- which.max(aggregate(values ~ nobs,
                                   data=plotData, mean)$values)
    print(paste("The observation with greater mean probability is:",
                maxProb))

    g <- ggplot(subset(plotData), aes(y=values, x=nobs)) + theme_bw()
    if (length(taus) > 1){
      if (scales.free) g <- g + facet_wrap(~taus, scales='free')
      else g <- g + facet_wrap(~taus)
    }

    g <- g + geom_point() +
      ylab("Posterior probability of being an outlier") + xlab("# Observation")

    if (plotProb) print(g)
  }

  return(plotData)
}
