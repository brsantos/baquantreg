#' Probability of being an outlier given the Bayesian quantile regression fits.
#'
#' Returns the probability that each observation is an outlier given the
#'  Bayesian quantile regression fits.
#'
#' @param object This is an object of the class "bqr", produced by a call
#'  to the bayesQR function.
#' @param burnin Initial part of the chain, which is to be discarded. Default
#'  value is 50.
#' @param plotKL If TRUE, the function prints the plot with all probabilities.
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

latentKL <- function(object, burnin = 50, plotKL = TRUE,
                     scales.free = FALSE, all.obs = TRUE, obs = 1) {

  if (class(object) != "bqr")
    stop("This function is not suited for your model.")

  taus <- object$tau
  nobs <- dim(object$chains[[1]]$vSample)[2]

  if (all.obs) seqObs <- 1:nobs
  else seqObs <- obs

  klValues <- sapply(object$chains, function(a){
    sapply(seqObs, function(b){
      otherV <- a$vSample[-c(1:burnin),-b]
      vSample <- a$vSample[-c(1:burnin), b]
      mean(sapply(1:(nobs-1), function(c){
        minV <- min(vSample, otherV[,c])
        maxV <- max(vSample, otherV[,c])
        g1 <- density(otherV[,c], from = minV, to = maxV)$y
        g2 <- density(vSample, from = minV, to = maxV)$y
        g1[g1 == 0] <- .Machine$double.eps
        g2[g2 == 0] <- .Machine$double.eps
        valF <- g2 * (log(g2) - log(g1))
        tail(cumsum(.5 * (valF[-1] + valF[-length(valF)])), 1)
      }))
    })
  })

  plotData <- data.frame(nobs = rep(seqObs, times=length(taus)),
                         values = as.numeric(klValues),
                         taus = rep(taus, each=length(seqObs)))

  if (all.obs){
    maxKL <- which.max(aggregate(values ~ nobs, data=plotData, mean)$values)

    print(paste("The observation with greater Kullback-Leibler divergence from the others is:",
                maxKL))

    g <- ggplot(subset(plotData), aes(y=values, x=nobs)) + theme_bw()
    if (length(taus) > 1){
      if (scales.free) g <- g + facet_wrap(~taus, scales='free')
      else g <- g + facet_wrap(~taus)
    }

    g <- g + geom_point() +
      ylab("Mean Kullback-Leibler divergence for each point") +
      xlab("# Observation")

    if (plotKL) print(g)
  }

  return(invisible(plotData))
}
