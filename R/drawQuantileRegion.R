#' Multiple-output Bayesian quantile regression model
#'
#' This function draws plots of the quantile region based on multiple-output
#'  quantile regression models.
#'
#' @param model This is an object of the class \code{multBQR}, produced by a
#'  call to the \code{multBayesQR} function.
#' @param datafile A data.frame from which to find the variables defined in the
#'  formula.
#' @param response Names of response variables
#' @param ngridpoints Number of grid points considered to build this quantile
#'  region, where a thorough search will look for the specified region, given
#'  the estimates for several directions. Default is 100, that will produce a
#'  grid with 10.000 points in the observed range of the data.
#' @param xValue Fixed value of the predictor variables. Default value is 1
#'  when there is only the intercept. If there is the interest in comparing the
#'  quantile regions for different values of predictors, it must used with a
#'  list of values for the predictors.
#' @param paintedArea If TRUE, it will plot the data points and the quantile
#'  region layer over the points, for each tau, in different plots. If FALSE,
#'  it will produce one plot showing showing all quantile regions for the
#'  different quantiles.
#' @param comparison Only considered when \code{paintedArea = FALSE}, if TRUE
#'  then it will plot comparisons of quantile regions for different values of
#'  \code{xvalue}.
#' @param result_folder Logical value determining whether all estimates are
#'  stored in a folder produced by a BayesX call. Default is FALSE.
#' @param path_folder If \code{result_folder = TRUE}, then one must inform the
#'  path where all results are stored.
#' @param splines_part  Logical value to indicate whether there are splines
#'  terms in the equation to draw the quantile contours.
#' @param wValue Fixed value to be plugged in the spline part of the equation.
#' @param print_plot Logical determining whether plot should be printed or
#'  data with coordinantes should be returned. Only checked when paintedArea is
#'  FALSE. Default is TRUE.
#' @param model_name When results will be collected in a folder, this should be
#'  the name of the name considered by BayesX to save all tables. Default is
#'  'bayesx.estim'.
#' @param name_var When there is a nonlinear variable from which one wants to
#'  consider different values for plotting, this should have the name of the
#'  variable.
#' @param ... Other parameters for \code{summary.multBQR}.
#' @return A ggplot with the quantile regions based on Bayesian quantile
#'  regression model estimates.
#' @useDynLib baquantreg

drawQuantileRegion <- function(model, datafile, response,
                               ngridpoints = 100, xValue = 1,
                               paintedArea = FALSE, comparison = FALSE,
                               result_folder = FALSE, path_folder = NULL,
                               splines_part = FALSE, wValue = NULL,
                               print_plot = TRUE,
                               model_name = 'bayesx.estim',
                               name_var, ...){

  if (!result_folder){
    directions <- sapply(model$modelsDir, function(a) a$direction)
    orthBases <- sapply(model$modelsDir, function(a) a$orthBasis)

    taus <- model$modelsDir[[1]]$tau
    ntaus <- length(taus)

    Y <- model$modelsDir[[1]]$data[ , model$response]

    estimates <- summary.multBQR(model, ...)

    if (model$method != 'bayesx'){
      betaEstimates <- lapply(estimates, function(a){
        sapply(a$BetaPosterior, function(b) b[ , 2])
      })

      betaDifDirections <- lapply(1:ntaus, function(a){
        sapply(betaEstimates, function(b) b[, a])
      })
    }
    else{
      betaEstimates <- lapply(estimates, function(a){
        sapply(a, function(b) b$BetaPosterior[, 1])
      })

      betaDifDirections <- lapply(1:ntaus, function(a){
        sapply(betaEstimates, function(b) b[, a])
      })
    }
  }
  else{
    if (is.null(path_folder))
      stop("You must define a path with all the results")
    else{
      results <- get_results(path_folder,
                             model_name = model_name,
                             splines = splines_part,
                             name_var = name_var)

      taus <- results$taus
      ntaus <- length(taus)

      Y <- datafile[, response]
      betaDifDirections <- results$betaDifDirections
      directions <- results$directions
      orthBases <- results$orthBases

      number_directions <- dim(directions)[2]

      splines_estimates <- results$spline_estimates_DifDirections
    }
  }

  y1range <- range(Y[,1])
  y2range <- range(Y[,2])

  seqY1 <- seq(y1range[1], y1range[2], length.out = ngridpoints)
  seqY2 <- seq(y2range[1], y2range[2], length.out = ngridpoints)

  pointsPlot <-  lapply(1:ntaus, function(a){
    if (!comparison){
      if (splines_part) spline_values <- sapply(1:number_directions,
                                                function(aa){
            estimates_direction <- splines_estimates[[a]][[aa]]
            distances <- abs(wValue - estimates_direction[, name_var])
            estimates_direction$pmean[which(distances == min(distances))[1]]
      })
      else spline_values <- rep(0, number_directions)

      checkPoints_values <- checkPoints(seqY1, seqY2, t(directions),
                                        t(orthBases), betaDifDirections[[a]],
                                        xValue, splines_part, spline_values)

      y1_ind <- which(rowSums(checkPoints_values) > 0)
      if(length(y1_ind) > 0){
        y1_inside <- seqY1[y1_ind]

        sub_checkPoints <- checkPoints_values[y1_ind, ]

        y2_inside_max <- apply(sub_checkPoints, 1, function(a)
          seqY2[max(which(a == 1))])
        y2_inside_min <- apply(sub_checkPoints, 1, function(a)
          seqY2[min(which(a == 1))])
      }
      else{
        y1_inside <- NA
        y2_inside_min <- NA
        y2_inside_max <- NA
      }

      data.frame(y1 = y1_inside,
                 min = y2_inside_min,
                 max = y2_inside_max)
    }
    else{
      lapply(xValue, function(b){
        checkPoints_values <- checkPoints(seqY1, seqY2,
                                          t(directions), t(orthBases),
                                          betaDifDirections[[a]],
                                          b)

        y1_ind <- which(rowSums(checkPoints_values) > 0)

        if(length(y1_ind) > 0){
          y1_inside <- seqY1[y1_ind]

          sub_checkPoints <- checkPoints_values[y1_ind, ]

          y2_inside_max <- apply(sub_checkPoints, 1, function(a)
            seqY2[max(which(a == 1))])
          y2_inside_min <- apply(sub_checkPoints, 1, function(a)
            seqY2[min(which(a == 1))])
        }
        else{
          y1_inside <- NA
          y2_inside_min <- NA
          y2_inside_max <- NA
        }

        data.frame(y1 = y1_inside,
                   min = y2_inside_min,
                   max = y2_inside_max)
      })
    }
  })

  if (paintedArea){
    g <- ggplot(data.frame(y1 = Y[,1], y2 = Y[,2])) + theme_bw()
    lapply(1:ntaus, function(a){
      g <- g + geom_point(aes(x = y1, y = y2))
      g + geom_ribbon(data = pointsPlot[[a]], aes(x = y1, ymin = min,
                                                  ymax = max),
                      alpha = 0.75, fill = 'lightblue') +
        ggtitle(paste0("Tau = ", taus[a]))
    })
  }
  else{
    if (comparison == FALSE){
      dataPlot <- data.frame(x = as.numeric(unlist(sapply(pointsPlot,
                                                          function(a){
        c(a$y1, rev(a$y1), a$y1[1])}))),
        y = as.numeric(unlist(sapply(pointsPlot, function(a){
          c(a$min, rev(a$max), a$min[1])}))),
        taus = as.numeric(unlist(sapply(1:ntaus, function(a){
          rep(taus[a],  2*dim(pointsPlot[[a]])[1] + 1)
        }))))

      if(print_plot){
        g <- ggplot() + theme_bw()
        g + geom_path(data = dataPlot, aes(x, y, linetype = factor(taus))) +
          scale_linetype_discrete(name = expression(tau)) +
          xlab(colnames(Y)[1]) +
          ylab(colnames(Y)[2]) +
          theme(legend.position = 'none')
      }
      else dataPlot
    }
    else{
      xPoints <- as.numeric(unlist(sapply(pointsPlot, sapply, function(a){
        c(a$y1, rev(a$y1), a$y1[1])
      })))

      yPoints <- as.numeric(unlist(sapply(pointsPlot, sapply, function(a){
        c(a$min, rev(a$max), a$min[1])})))

      tausPoints <- as.numeric(unlist(lapply(1:ntaus, function(a){
        sapply(1:length(xValue), function(b){
          rep(taus[a], 2*dim(pointsPlot[[a]][[b]])[1] + 1)
      })})))

      predictorsType <- as.numeric(unlist(lapply(1:ntaus, function(a){
        sapply(1:length(xValue), function(b){
          rep(b, 2*dim(pointsPlot[[a]][[b]])[1] + 1)
      })})))

      dataPlot <- data.frame(x = xPoints,
                             y = yPoints,
                             taus = tausPoints,
                             predictors = predictorsType)

      if(print_plot){
        g <- ggplot() + theme_bw()
        g + geom_path(data = dataPlot, aes(x, y, linetype = factor(taus),
                                           color = factor(predictors))) +
          scale_linetype_discrete(name = expression(tau)) +
          xlab(colnames(Y)[1]) +
          ylab(colnames(Y)[2]) +
          theme(legend.position = 'none')
      }
      else dataPlot
    }
  }
}