#' Multiple-output Bayesian quantile regression model
#'
#' This function checks whether points belong to quantile regions or not, based
#'  on estimated models for models with 2 dimensions.
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
#' @param lambda_a Parameter to be make the adjustment of the intercept to
#'  correct the miscoverage of the quantile regions.
#' @param ... Other parameters for \code{summary.multBQR}.
#' @return A ggplot with the quantile regions based on Bayesian quantile
#'  regression model estimates.
#' @useDynLib baquantreg

checkpoints_qreg_2D <- function(model, datafile, response,
                                points_y, x_values = 1,
                                path_folder = NULL,
                                splines_part = FALSE, w_values = NULL,
                                model_name = 'bayesx.estim',
                                name_var, lambda_a, ...){

  if (is.null(path_folder))
    stop("You must define a path with all the results")
  else {
    results <- get_results(path_folder,
                           model_name = model_name,
                           splines = splines_part,
                           name_var = name_var,
                           n_dim = 2, datafile, response,
                           lambda_a = lambda_a)

    taus <- results$taus
    ntaus <- length(taus)

    Y <- datafile[, response]

    betaDifDirections <- results$betaDifDirections
    splines_estimates <- results$spline_estimates_DifDirections

    # if (!upperq){
    #   betaDifDirections <- results$betaDifDirections
    #   splines_estimates <- results$spline_estimates_DifDirections
    # }
    # else {
    #   betaDifDirections <- results$upperq_DifDirections
    #   splines_estimates <- results$upperq_spline_estimates
    # }
    #
    # if (!lowerq){
    #   betaDifDirections <- results$betaDifDirections
    #   splines_estimates <- results$spline_estimates_DifDirections
    # }
    # else {
    #   betaDifDirections <- results$lowerq_DifDirections
    #   splines_estimates <- results$lowerq_spline_estimates
    # }

    directions <- results$directions
    orthBases <- results$orthBases

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
                                        t(directions),
                                        t(orthBases),
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