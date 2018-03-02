#' Plot methods for summary of Bayesian quantile regression models.
#'
#' Returns a ggplot of the estimates of the models and their respective
#'  credible intervals.
#'
#' @param x This is an object of the class "summary.zitobitQR",
#' produced by a call to the summary.zitobitQR function.
#' @param separate if FALSE, all plots are returned in the same output, and
#'  if TRUE the plot for each variable is returned separately. Default value
#'  is FALSE.
#' @param beta If TRUE, it plots the posterior estimates for beta(tau).
#'  Default is set to TRUE.
#' @param sigma If TRUE, it plots the posterior estimates for sigma.
#'  Default is set to FALSE.
#' @param gamma If TRUE, it plots the posterior estimates for gamma.
#'  Default is set to FALSE.
#' @param ... other plot params (currently not used)
#' @return A ggplot of the posterior estimates with their credible intervals.
#' @export
#' @useDynLib baquantreg
#' @import ggplot2

plot.summary.zitobitQR <- function(x, separate = FALSE, beta = TRUE,
                                   sigma = FALSE, gamma = FALSE, ...){
  if (class(x) != "summary.zitobitQR")
    stop("Look for the correct method for your model.")

  ntaus <- length(x$BetaPosterior)
  nvar <- dim(x$BetaPosterior[[1]])[1]

  vnames <- x$BetaPosterior[[1]][,1]

  postEst <- as.numeric(sapply(x$BetaPosterior, function(a) a[,2]))
  lowerQ <- as.numeric(sapply(x$BetaPosterior, function(a) a[,3]))
  upperQ <- as.numeric(sapply(x$BetaPosterior, function(a) a[,4]))

  postEst_g <- as.numeric(sapply(x$GammaPosterior, function(a) a[,2]))
  lowerQ_g <- as.numeric(sapply(x$GammaPosterior, function(a) a[,3]))
  upperQ_g <- as.numeric(sapply(x$GammaPosterior, function(a) a[,4]))

  taus <- x$taus

  plotData <- data.frame(vnames = rep(vnames, times=ntaus),
                         taus = rep(taus, each=nvar),
                         postEst = postEst,
                         lowerQ = lowerQ,
                         upperQ = upperQ)

  plotDataGamma <- data.frame(vnames = rep(vnames, times=ntaus),
                              taus = rep(taus, each=nvar),
                              postEst = postEst_g,
                              lowerQ = lowerQ_g,
                              upperQ = upperQ_g)

  if (sigma){
    g1 <- ggplot(x$SigmaPosterior, aes(x=taus)) + theme_bw()
    g1 <- g1 + geom_line(aes(y=Mean)) +
      geom_line(aes(y=Lower), linetype=2) +
      geom_line(aes(y=Upper), linetype=2) +
      ylab("Posterior estimates for sigma") +
      xlab(expression(tau))
  }
  if (sigma) print(g1)

  if (beta){
    if (!separate){
      g <- ggplot(plotData, aes(x=taus)) + theme_bw()
      g <- g + geom_line(aes(y=postEst)) + geom_line(aes(y=lowerQ), linetype=2) +
        ylab("Posterior estimates") + xlab(expression(tau)) +
        geom_line(aes(y=upperQ), linetype=2) + facet_wrap(~vnames, scales='free')

    }
    else {
      lapply(vnames, function(a){
        g <- ggplot(subset(plotData, vnames==a), aes(x=taus)) + theme_bw()
        g <- g + geom_line(aes(y=postEst)) + geom_line(aes(y=lowerQ), linetype=2) +
          ylab("Posterior estimates") + xlab(expression(tau)) +
          geom_line(aes(y=upperQ), linetype=2) + facet_wrap(~vnames,
                                                            scales='free')
      })
    }
  }
  if (beta) print(g)

  if (gamma){
    if (!separate){
      gG <- ggplot(plotDataGamma, aes(x=taus)) + theme_bw()
      gG <- gG + geom_line(aes(y=postEst)) + geom_line(aes(y=lowerQ), linetype=2) +
        ylab("Posterior estimates") + xlab(expression(tau)) +
        geom_line(aes(y=upperQ), linetype=2) + facet_wrap(~vnames, scales='free')

    }
    else {
      lapply(vnames, function(a){
        gG <- ggplot(subset(plotDataGamma, vnames==a), aes(x=taus)) + theme_bw()
        gG <- gG +
          geom_line(aes(y=postEst))+geom_line(aes(y=lowerQ), linetype=2) +
          ylab("Posterior estimates") + xlab(expression(tau)) +
          geom_line(aes(y=upperQ), linetype=2) +
          facet_wrap(~vnames, scales='free')
      })
    }
  }
  if (gamma) print(gG)
}