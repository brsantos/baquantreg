#' Plot function to compare the posterior distribution of the censoring
#'  indicator given a factor used in the model.
#'
#' Returns a plot comparing the posterior distribution for
#'
#' @param object This is an object of the class "zitobitQR", produced by a call
#'  to the zitobitQR function.
#' @param variable Name of the variable to be used in the comparison.
#' @return A plot of the posterior density of the probability of being censored
#'  given some factor used in the model.
#' @export
#' @useDynLib baquantreg
#' @import ggplot2

plotCensoring <- function(object, variable){
  if (class(object) != "zitobitQR")
    stop("This function is not suited for your model.")

  taus <- object$tau

  y <- object$y
  varInd <- as.numeric(y == 0)
  valuesProb <- as.numeric(sapply(object$chains, function(a){
    a$indCens
  }))

  varPlot <- as.factor(object$data[, variable])

  plotData <- data.frame(values = valuesProb,
                         taus = rep(taus, length(y)),
                         varPlot = rep(varPlot, times = length(taus)),
                         varInd = rep(varInd, times = length(taus)))

  plotData <- subset(plotData, varInd == 1,
                     select=c("values", "taus", "varPlot"))

  g <- ggplot(subset(plotData), aes(y=..density.., x=values)) + theme_bw()
  if (length(taus) > 1) g <- g + facet_wrap(~taus, scales='free')

  g <- g + geom_density(aes(linetype=varPlot)) +
    ylab("Posterior distribution for the censoring probability") +
    xlab("") + theme(legend.position = 'bottom') +
    scale_linetype(name = variable)

  print(g)
  return(invisible(plotData))
}