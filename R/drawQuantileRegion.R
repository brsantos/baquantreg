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
#' @param lower_q Logical determining whether quantile region should be based
#'  on the lower value of 0.95 credible interval for each parameter of the model.
#' @param upper_q Logical determining whether quantile region should be based
#'  on the upper value of 0.95 credible interval for each parameter of the model.
#' @param lambda_a Parameter to be make the adjustment of the intercept to
#'  correct the miscoverage of the quantile regions.
#'  \code{lower_q} must be FALSE for this argument to be checked.
#' @param range_y matrix type object containing in the first line the range for
#'  the first dimension and in the second line the range for the second
#'  dimension of Y. This will be used to find the respective quantile regions.
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
                               name_var, range_y = NULL, lower_q = FALSE, upper_q = FALSE,
                               lambda_a = NULL,
                               ...){

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
    else {
      results <- get_results(path_folder,
                             model_name = model_name,
                             splines = splines_part,
                             name_var = name_var,
                             n_dim = 2,
                             datafile,
                             response,
                             lambda_a = lambda_a)

      taus <- results$taus
      ntaus <- length(taus)

      Y <- datafile[, response]
      directions <- results$directions
      orthBases <- results$orthBases

      number_directions <- dim(directions)[2]

      if (!any(c(lower_q, upper_q))){
        betaDifDirections <- results$betaDifDirections
        splines_estimates <- results$spline_estimates_DifDirections
      }
      else if (lower_q){
        betaDifDirections <- lapply(1:length(results$betaDifDirections),
                                    function(xxx){
          results$betaDifDirections[[xxx]] -
            stats::qnorm(0.975) * diag(c(xValue, 0)) %*%
            results$sdDifDirections[[xxx]]
        })
        splines_estimates <- results$spline_estimates_DifDirections
      }
      else {
        betaDifDirections <- lapply(1:length(results$betaDifDirections),
                                    function(xxx){
          results$betaDifDirections[[xxx]] +
            stats::qnorm(0.975) * diag(c(xValue, 0)) %*%
            results$sdDifDirections[[xxx]]
        })
        splines_estimates <- results$spline_estimates_DifDirections
      }
    }
  }

  if(is.null(range_y)){
    y1range <- range(Y[,1])
    y2range <- range(Y[,2])
  }
  else {
    y1range <- range_y[1, ]
    y2range <- range_y[2, ]
  }

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

      if(length(checkPoints_values) > 0){
        y1_inside <- unique(sort(checkPoints_values[, 1]))

        y2_inside_max <- sapply(y1_inside, function(a){
          values_to_filter <- checkPoints_values[ , 1] == a
          max(checkPoints_values[ values_to_filter, 2])
        })

        y2_inside_min <- sapply(y1_inside, function(a){
          values_to_filter <- checkPoints_values[ , 1] == a
          min(checkPoints_values[values_to_filter, 2])
        })
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

        if(length(checkPoints_values) > 0){
          y1_inside <- unique(sort(checkPoints_values[, 1]))

          y2_inside_max <- sapply(y1_inside, function(a){
            values_to_filter <- checkPoints_values[ , 1] == a
            max(checkPoints_values[ values_to_filter, 2])
          })

          y2_inside_min <- sapply(y1_inside, function(a){
            values_to_filter <- checkPoints_values[ , 1] == a
            min(checkPoints_values[values_to_filter, 2])
          })
        }
        else {
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

      if (print_plot){
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