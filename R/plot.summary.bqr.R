#' Plot methods for summary of Bayesian quantile regression models.
#'
#' Returns a ggplot of the estimates of the models and their respective
#'  credible intervals.
#'
#' @param x This is an object of the class \code{summary.bqr},
#' produced by a call to the summary.bqr function.
#' @param separate if FALSE, all plots are returned in the same output, and
#'  if TRUE the plot for each variable is returned separately. Default value
#'  is FALSE.
#' @param sigma if FALSE, it does the plot the posterior estimates for quantile
#'  regression parameters, and if TRUE, it plots the posterior estimates only
#'  for sigma. Default is set to FALSE.
#' @param ... Other \code{plot} arguments.
#' @return A ggplot of the posterior estimates with their credible intervals.
#' @export
#' @useDynLib baquantreg
#' @import ggplot2

plot.summary.bqr <- function(x, separate = F, sigma = F, ...){
  if (class(x) != "summary.bqr")
    stop("Look for the correct method for your model.")

  ntaus <- length(x$BetaPosterior)
  nvar <- dim(x$BetaPosterior[[1]])[1]

  vnames <- x$BetaPosterior[[1]][,1]

  postEst <- as.numeric(sapply(x$BetaPosterior, function(a) a[,2]))
  lowerQ <- as.numeric(sapply(x$BetaPosterior, function(a) a[,3]))
  upperQ <- as.numeric(sapply(x$BetaPosterior, function(a) a[,4]))

  taus <- x$taus

  plotData <- data.frame(vnames = rep(vnames, times=ntaus),
                         taus = rep(taus, each=nvar),
                         postEst = postEst,
                         lowerQ = lowerQ,
                         upperQ = upperQ)

  if (sigma){
    g1 <- ggplot(x$SigmaPosterior, aes(x=taus)) + theme_bw()
    g1 + geom_line(aes(y=Mean)) +
      geom_line(aes(y=Lower), linetype=2) +
      geom_line(aes(y=Upper), linetype=2) +
      ylab("Posterior estimates for sigma") +
      xlab(expression(tau))
  }
  else {
    if (!separate){
      g <- ggplot(plotData, aes(x=taus)) + theme_bw()
      g + geom_line(aes(y=postEst)) + geom_line(aes(y=lowerQ), linetype=2) +
        ylab("Posterior estimates") + xlab(expression(tau)) +
        geom_line(aes(y=upperQ), linetype=2) + facet_wrap(~vnames, scales='free')

    }
    else {
      lapply(vnames, function(a){
        g <- ggplot(subset(plotData, vnames==a), aes(x=taus)) + theme_bw()
        g + geom_line(aes(y=postEst)) + geom_line(aes(y=lowerQ), linetype=2) +
          ylab("Posterior estimates") + xlab(expression(tau)) +
          geom_line(aes(y=upperQ), linetype=2) + facet_wrap(~vnames,
                                                            scales='free')
      })
    }
  }
}