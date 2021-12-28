#' Multiple-output Bayesian quantile regression model
#'
#' This function checks whether points belong to quantile regions or not, based
#'  on estimated models for models with 4 dimensions.
#'
#' @param model This is an object of the class \code{multBQR}, produced by a
#'  call to the \code{multBayesQR} function.
#' @param datafile A data.frame from which to find the variables defined in the
#'  formula.
#' @param response Names of response variables
#' @param points_y the exact points in Y, in which one wants to find its
#'  respective quantile region.
#' @param x_values Fixed value of the predictor variables.
#' @param path_folder The path where all results are stored.
#' @param splines_part  Logical value to indicate whether there are splines
#'  terms in the equation to draw the quantile contours.
#' @param w_values Value to be considered in the nonlinear part of the model.
#' @param model_name When results will be collected in a folder, this should be
#'  the name of the name considered by BayesX to save all tables. Default is
#'  'bayesx.estim'.
#' @param name_var When there is a nonlinear variable from which one wants to
#'  consider different values for plotting, this should have the name of the
#'  variable.
#' @param adaptive_dir If  \code{TRUE}, then directions will take into account
#'  the marginal quantiles of each dimension of the response variable.
#'  Otherwise, the direction vector are created creating all possible
#'  combinations of points inside the interval [-1, 1] given the number of
#'  points \code{directionPoint}. The default is \code{FALSE}.
#' @param lowerq If TRUE, then it will take into account the lower quantiles
#'  for each coefficient of the model.
#' @param upperq If TRUE, then it will take into account the upper quantiles
#'  for each coefficient of the model.
#' @param ... Other parameters for \code{summary.multBQR}.
#' @return A ggplot with the quantile regions based on Bayesian quantile
#'  regression model estimates.
#' @useDynLib baquantreg

checkpoints_qreg_4D <- function(model, datafile, response,
                                points_y, x_values = 1,
                                path_folder = NULL,
                                splines_part = FALSE, w_values = NULL,
                                model_name = 'bayesx.estim',
                                name_var, adaptive_dir = FALSE,
                                upperq = FALSE,
                                lowerq = FALSE, ...){

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
    if (!upperq){
      betaDifDirections <- results$betaDifDirections
      splines_estimates <- results$spline_estimates_DifDirections
    }
    else {
      betaDifDirections <- results$upperq_DifDirections
      splines_estimates <- results$upperq_spline_estimates
    }

    if (!lowerq){
      betaDifDirections <- results$betaDifDirections
      splines_estimates <- results$spline_estimates_DifDirections
    }
    else {
      betaDifDirections <- results$lowerq_DifDirections
      splines_estimates <- results$lowerq_spline_estimates
    }

    directions <- results$directions
    orthBases1 <- results$orthBases1
    orthBases2 <- results$orthBases2
    orthBases3 <- results$orthBases3

    number_directions <- dim(directions)[2]
  }

  n_points <- dim(points_y)[1]

  points_inside <-  lapply(1:ntaus, function(a){
    sapply(1:n_points, function(nnn){
      if (splines_part){
        spline_values <-
          sapply(1:number_directions, function(aa){
            estimates_direction <- splines_estimates[[a]][[aa]]
            distances <- abs(w_values[nnn] - estimates_direction[, name_var])
            estimates_direction[, 2][which(distances == min(distances))[1]]
          })
      } else spline_values <- rep(0, number_directions)

      checkPoints_val <- checkPoints_4d(points_y[nnn, 1],
                                        points_y[nnn, 2],
                                        points_y[nnn, 3],
                                        points_y[nnn, 4],
                                        t(directions),
                                        t(orthBases1),
                                        t(orthBases2),
                                        t(orthBases3),
                                        betaDifDirections[[a]],
                                        x_values[nnn, ],
                                        splines_part, spline_values)
      if(length(checkPoints_val) == 0) inside <- FALSE
      else inside <- TRUE

      inside
    })
  })
  points_inside
}