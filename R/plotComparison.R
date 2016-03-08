#' Plot function to compare the posterior distribution for the latent variables
#'
#' Returns a plot comparing the posterior distribution for an observation and
#'  the distribution for all other observations.
#'
#' @param object This is an object of the class "bqr", produced by a call
#'  to the bayesQR function.
#' @param observation Number of the observation to be compared with all the
#'  other observations, in the posterior distribution of their latent
#'  variables.
#' @param burnin Initial part of the chain, which is to be discarded. Default
#'  value is 50.
#' @param plotComp If TRUE, the function prints the plot with the comparison.
#'  Default is set to TRUE.
#' @return A plot of the posterior density of the latent variable of a selected
#'  observation and the correspondent density for all other observations.
#' @export
#' @useDynLib baquantreg
#' @import ggplot2

plotComparison <- function(object, observation, burnin = 50, plotComp = T){
  if (class(object) != "bqr")
    stop("This function is not suited for your model.")

  taus <- object$tau

  dataAll <- as.numeric(sapply(object$chains, function(a){
    apply(a$vSample[-c(1:burnin),-observation], 2, mean)
  }))

  dataObs <- as.numeric(sapply(object$chains, function(a){
    as.numeric(a$vSample[-c(1:burnin),observation])
  }))

  t1 <- length(dataAll)/length(taus)
  t2 <- length(dataObs)/length(taus)

  plotData <- data.frame(values = c(dataAll, dataObs),
                         taus = c(rep(taus, each=t1),
                                  rep(taus, each=t2)),
                         type = c(rep("Others", times=length(dataAll)),
                                  rep("Observation", times=length(dataObs))))

  g <- ggplot(subset(plotData), aes(y=..density..)) + theme_bw()
  if (length(taus) > 1) g <- g + facet_wrap(~taus, scales='free')

  g <- g + geom_density(aes(linetype=type, x=values)) +
    ylab("Posterior distribution for the latent variables") +
    xlab("") + theme(legend.position = 'bottom')

  if (plotComp) print(g)
  return(plotData)
}