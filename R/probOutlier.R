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

probOutlier <- function(object, burnin = 50, plotProb = TRUE,
                        scales.free = FALSE){

  if (class(object) != "bqr")
    stop("This function is not suited for your model.")

  taus <- object$tau
  nobs <- dim(object$chains[[1]]$vSample)[2]

  prob <- t(sapply(1:nobs, function(b){
    dataAll <- as.numeric(sapply(object$chains, function(a){
      apply(a$vSample[-c(1:burnin),-b], 2, mean)
    }))

    dataObs <- as.numeric(sapply(object$chains, function(a){
      as.numeric(a$vSample[-c(1:burnin),b])
    }))

    t1 <- length(dataAll)/length(taus)
    t2 <- length(dataObs)/length(taus)

    probData <- data.frame(values = c(dataAll, dataObs),
                           taus = c(rep(taus, each=t1),
                                    rep(taus, each=t2)),
                           type = c(rep("Others", times=length(dataAll)),
                                    rep("Observation", times=length(dataObs))))

    sapply(taus, function(c){
      sum( subset(probData, taus==c & type == "Observation")$values >
            max(subset(probData, taus==c & type != "Observation")$values) ) /
        length(subset(probData, taus==c & type == "Observation")$values)
    })
  }))

  plotData <- data.frame(nobs = rep(1:nobs, times=dim(prob)[2]),
                         values = as.numeric(prob),
                         taus = rep(taus, each=dim(prob)[1]))

  maxProb <- which.max(aggregate(values ~ nobs, data=plotData, mean)$values)
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

  return(plotData)
}
