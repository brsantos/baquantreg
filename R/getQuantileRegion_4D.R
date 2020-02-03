#' Multiple-output Bayesian quantile regression model
#'
#' This function draws plots of the quantile region based on multiple-output
#'  quantile regression models, when the response variable has 4 dimensions.
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
#' @param path_folder The path where all results are stored.
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
#' @param adaptive_dir If  \code{TRUE}, then directions will take into account the
#'  marginal quantiles of each dimension of the response variable. Otherwise,
#'  the direction vector are created creating all possible combinations of
#'  points inside the interval [-1, 1] given the number of points
#'  \code{directionPoint}. The default is \code{FALSE}.
#' @param ... Other parameters for \code{summary.multBQR}.
#' @return A ggplot with the quantile regions based on Bayesian quantile
#'  regression model estimates.
#' @useDynLib baquantreg

getQuantileRegion_4D <- function(model, datafile, response,
                                  ngridpoints = 100, xValue = 1,
                                  path_folder = NULL,
                                  splines_part = FALSE, wValue = NULL,
                                  print_plot = TRUE,
                                  model_name = 'bayesx.estim',
                                  name_var, adaptive_dir = FALSE, ...){

  if (is.null(path_folder))
    stop("You must define a path with all the results")
  else {
    results <- get_results(path_folder,
                           model_name = model_name,
                           splines = splines_part,
                           name_var = name_var,
                           n_dim = 4, datafile, response, adaptive_dir)

    taus <- results$taus
    ntaus <- length(taus)

    Y <- datafile[, response]
    betaDifDirections <- results$betaDifDirections
    directions <- results$directions
    orthBases1 <- results$orthBases1
    orthBases2 <- results$orthBases2
    orthBases3 <- results$orthBases3

    number_directions <- dim(directions)[2]

    splines_estimates <- results$spline_estimates_DifDirections
  }

  y1range <- range(Y[, 1])
  y2range <- range(Y[, 2])
  y3range <- range(Y[, 3])
  y4range <- range(Y[, 4])

  seqY1 <- seq(y1range[1], y1range[2], length.out = ngridpoints)
  seqY2 <- seq(y2range[1], y2range[2], length.out = ngridpoints)
  seqY3 <- seq(y3range[1], y3range[2], length.out = ngridpoints)
  seqY4 <- seq(y4range[1], y4range[2], length.out = ngridpoints)

  pointsPlot <-  lapply(1:ntaus, function(a){
    if (splines_part){
      spline_values <-
        sapply(1:number_directions, function(aa){
          estimates_direction <- splines_estimates[[a]][[aa]]
          distances <- abs(wValue - estimates_direction[, name_var])
          estimates_direction$pmean[which(distances == min(distances))[1]]
        })
    } else spline_values <- rep(0, number_directions)

    checkPoints_val <- checkPoints_4d(seqY1, seqY2, seqY3, seqY4,
                                         t(directions),
                                         t(orthBases1), t(orthBases2),
                                         t(orthBases3),
                                         betaDifDirections[[a]],
                                         xValue, splines_part, spline_values)

    all_points <- lapply(unique(sort(checkPoints_val[, 4])), function(dddd){
      checkPoints_values <- checkPoints_val[checkPoints_val[,4] == dddd, ]

      y3_inside <- unique(checkPoints_values[, 3])

      points_inside <- lapply(y3_inside, function(aaa){
        checkPoints_values_aux_ind <- checkPoints_values[,3] == aaa
        checkPoints_values_aux <-
          checkPoints_values[checkPoints_values_aux_ind, 1:2]

        if(length(checkPoints_values_aux) > 0){
          y1_inside <- unique(sort(checkPoints_values[, 1]))

          y2_inside_max <- sapply(y1_inside, function(a){
            values_to_filter <- checkPoints_values_aux[ , 1] == a
            max(checkPoints_values_aux[ values_to_filter, 2])
          })

          y2_inside_min <- sapply(y1_inside, function(a){
            values_to_filter <- checkPoints_values_aux[ , 1] == a
            min(checkPoints_values_aux[values_to_filter, 2])
          })

        }
        else {
          y1_inside <- NA
          y2_inside_min <- NA
          y2_inside_max <- NA
        }

        data.frame(y1 = rep(y1_inside, times = 2),
                   y2 = c(y2_inside_min, y2_inside_max),
                   y3 = rep(aaa, length(y1_inside) * 2),
                   type = rep(c('min', 'max'), each = length(y1_inside)))
      })

      points_inside_all <- do.call(rbind.data.frame, points_inside)
      points_inside_all[stats::complete.cases(points_inside_all), ]
      points_inside_all$y4 <- rep(dddd, dim(points_inside_all)[1])
    })
    all_points_all <- do.call(rbind.data.frame, all_points)
    all_points_all[stats::complete.cases(all_points_all), ]
  })
}