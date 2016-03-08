#' Plot function to analyse the posterior latent variables in a Bayesian
#' quantile regression model.
#'
#' Returns a ggplot of the posterior trace of the latent variables
#'
#' @param object This is an object of the class "bqr",
#' produced by a call to the bayesQR function.
#' @param separate if FALSE, all plots are returned in the same output, and
#'  if TRUE the plot for each variable is returned separately. Default value
#'  is FALSE.
#' @return A ggplot of the posterior trace for all latent variables.
#' @export
#' @useDynLib baquantreg
#' @import ggplot2

plotvfit <- function(object, separate = F, burnin = 0, ...){
  if (class(object) != "bqr")
    stop("This function is not suited for your model.")

  taus <- object$tau

  dataV <- object$chains[[1]]$vSample
  nobs <- dim(object$chains[[1]]$vSample)[2]
  sizeChain <- dim(object$chains[[1]]$vSample)[1]

  if (length(taus) > 1){
    for (i in 2:length(taus))
      dataV <- rbind(dataV, object$chains[[i]]$vSample)
  }

  plotData <- data.frame(dataV,
                         taus = rep(taus, each=sizeChain),
                         sizeChain = rep(1:sizeChain, times=length(taus)))

  if (!separate){
    g <- ggplot(subset(plotData, sizeChain > burnin),
                aes(x=sizeChain)) + theme_bw()
    if (length(taus) > 1) g <- g + facet_wrap(~taus, scales='free')

    for (i in 1:nobs){
      g <- g + geom_line(aes_string(y=paste("X", i, sep="")),
                         colour='grey50', alpha=0.1) +
        ylab("Trace for latent variables") +
        xlab("Iteration")
    }
    print(g)
  }
  else {
    lapply(taus, function(a){
      g <- ggplot(subset(plotData, taus==a & sizeChain > burnin),
                  aes(x=sizeChain)) + theme_bw()
      for (i in 1:nobs){
        g <- g + geom_line(aes_string(y=paste("X", i, sep="")),
                           colour='grey50', alpha=0.1) +
          ylab("Trace for latent variables") +
          xlab("Iteration") + facet_wrap(~taus, scales='free')
      }
      print(g)
    })
  }
}