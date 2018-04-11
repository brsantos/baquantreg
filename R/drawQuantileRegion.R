#' Multiple-output Bayesian quantile regression model
#'
#' This function draws plots of the quantile region based on multiple-output quantile
#' regression models.
#'
#' @param model This is an object of the class \code{multBQR}, produced by a call to the
#'  \code{multBayesQR} function.
#' @param ngridpoints Number of grid points considered to build this quantile region,
#' where a thorough search will look for the specified region, given the estimates
#' for several directions. Default is 100, that will produce a grid with 10.000 points
#' in the observed range of the data.
#' @param xValue Fixed value of the predictor variables. It is not necessary for models
#' with only the intercept.
#' @param paintedArea If TRUE, it will plot the data points and the quantile region layer
#' over the points, for each tau, in different plots. If FALSE, it will produce one
#' plot showing showing all quantile regions for the different quantiles.
#' @param ... Other parameters for \code{summary.multBQR}.
#' @return A ggplot with the quantile regions based on Bayesian quantile regression
#' model estimates.
#' @useDynLib baquantreg
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise


drawQuantileRegion <- function(model, ngridpoints = 100, xValue, paintedArea = TRUE, ...){
  estimates <- summary.multBQR(model, ...)

  taus <- model[[1]]$tau
  ntaus <- length(taus)

  Y <- stats::model.extract(stats::model.frame(model[[1]]$formula, model[[1]]$data), "response")

  directions <- sapply(model, function(a) a$direction)
  orthBases <- sapply(model, function(a) a$orthBasis)

  y1range <- range(Y[,1])
  y2range <- range(Y[,2])

  seqY1 <- seq(y1range[1], y1range[2], length.out = ngridpoints)
  seqY2 <- seq(y2range[1], y2range[2], length.out = ngridpoints)

  Yseq <- cbind(rep(seqY1, times = ngridpoints),
                rep(seqY2, each = ngridpoints))

  betaEstimates <- lapply(estimates, function(a){
    sapply(a$BetaPosterior, function(b) b[ , 2])
  })

  betaDifDirections <- lapply(1:ntaus, function(a){
    sapply(betaEstimates, function(b) b[,a])
  })

  YResp <- Yseq %*% directions

  Xdirection <- Yseq %*% orthBases

  pointsPlot <-  lapply(1:ntaus, function(a){
    if(dim(betaDifDirections[[a]])[1] == 2){
      linPredictor <- matrix(rep(betaDifDirections[[a]][1,], ngridpoints^2), ncol = ncol(directions), byrow = TRUE) +
        Xdirection * matrix(rep(betaDifDirections[[a]][2,], ngridpoints^2), ncol = ncol(directions), byrow = TRUE)
    }
    else{
      linPredictor <- matrix(rep(betaDifDirections[[a]][1,], ngridpoints^2), ncol = ncol(directions), byrow = TRUE) +
        matrix(rep(t(betaDifDirections[[a]][2:(dim(betaDifDirections[[a]])[1] - 1), ]) %*% xValue,
                   ngridpoints^2),
               ncol = ncol(directions), byrow = TRUE) +
        Xdirection * matrix(rep(betaDifDirections[[a]][dim(betaDifDirections[[a]])[1], ], ngridpoints^2),
                            ncol = ncol(directions), byrow = TRUE)
    }

    matIndices <- YResp > linPredictor
    indices <- apply(matIndices, 1, all)

    dataRegion <- data.frame(y1 = Yseq[indices, 1],
                             y2 = Yseq[indices, 2])

    summarize(group_by(dataRegion, y1), min = min(y2),
                     max = max(y2))

  })

  if (paintedArea){
    g <- ggplot(data.frame(y1 = Y[,1], y2 = Y[,2])) + theme_bw()
    lapply(1:ntaus, function(a){
      g <- g + geom_point(aes(x = y1, y = y2))
      g + geom_ribbon(data = pointsPlot[[a]], aes(x = y1, ymin = min, ymax = max),
                      alpha = 0.75, fill = 'lightblue') +
        ggtitle(paste0("Tau = ", taus[a]))
    })
  }
  else{
    dataPlot <- data.frame(x = as.numeric(unlist(sapply(pointsPlot, function(a){
      c(a$y1, rev(a$y1), a$y1[1])}))),
      y = as.numeric(unlist(sapply(pointsPlot, function(a){
        c(a$min, rev(a$max), a$min[1])}))),
      taus = as.numeric(unlist(sapply(1:ntaus, function(a){
        rep(taus[a],  2*dim(pointsPlot[[a]])[1] + 1)
      }))))

    g <- ggplot(dataPlot) + theme_bw()
    g + geom_path(aes(x, y, linetype = factor(taus))) +
      scale_linetype_discrete(name = expression(tau)) +
      xlab(colnames(Y)[1]) +
      ylab(colnames(Y)[2]) +
      theme(legend.position = 'none') + xlim(y1range) + ylim(y2range)
  }
}