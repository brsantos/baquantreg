## Organizing the results in a folder

get_results <- function(path_folder, model_name, splines = FALSE, name_var,
                        n_dim = 2, datafile_original, response, adaptive_dir,
                        lambda_a = NULL){
  folders <- list.files(path_folder)

  ## Considering at most 999 directions
  length_directions <- attr(regexpr("dir_[0-9]+_tau", folders), 'match.length')
  directions_ind <- as.numeric(ifelse(length_directions == 10,
                                      substr(folders, 5, 6),
                                      ifelse(length_directions == 9,
                                             substr(folders, 5, 5),
                                             substr(folders, 5, 7))))

  # Considering taus with at most two decimal places
  start_positions <- regexpr("tau_0.[0-9]+", folders)
  length_tau <- attr(regexpr("tau_0.[0-9]+", folders), 'match.length')
  taus <- as.numeric(ifelse(length_tau == 7, substr(folders, start_positions+4,
                                                    start_positions + 6),
                            substr(folders, start_positions + 4,
                                   start_positions + 7)))
  unique_taus <- unique(taus)

  ## Going through all folders and getting all estimates and draws.
  results <- lapply(folders, function(a){
    fixed_effects_files <- grepl("_FixedEffects[0-9]+.res",
                                 list.files(paste0(path_folder, '/', a, '/')))
    beta_draws_files <- grepl("_FixedEffects[0-9]+_sample.raw",
                              list.files(paste0(path_folder, '/', a, '/')))

    files_results <- list.files(paste0(path_folder, '/', a,
                                       '/'))[fixed_effects_files]
    draws_results <- list.files(paste0(path_folder, '/', a,
                                       '/'))[beta_draws_files]

    if(sum(fixed_effects_files) > 1){
      all_files <- lapply(files_results, function(aa){
        utils::read.table(paste0(path_folder, '/', a, '/', aa), head = TRUE)
      })
      fixedEffects <- do.call(rbind, all_files)[, 3]
      fixedEffects_sd <- do.call(rbind, all_files)[, 4]
      fixedEffects_lowerq <- do.call(rbind, all_files)[, 5]
      fixedEffects_upperq <- do.call(rbind, all_files)[, 9]
      varnames <- do.call(rbind, all_files)[, 2]

      if (!is.null(lambda_a)){
        fixedEffects[1] <- fixedEffects[1] - lambda_a
      }
    } else {
      info <- utils::read.table(paste0(path_folder, '/', a, '/',
                                       files_results), head = TRUE)
      fixedEffects <- info[, 3]
      fixedEffects_sd <- info[, 4]
      fixedEffects_lowerq <- info[, 5]
      fixedEffects_upperq <- info[, 9]
      varnames <- info[, 2]

      if (!is.null(lambda_a)){
        fixedEffects[1] <- fixedEffects[1] - lambda_a
      }
    }

    if (sum(beta_draws_files) > 1){
      all_files <- lapply(draws_results, function(aa){
        utils::read.table(paste0(path_folder, '/', a, '/', aa), head = TRUE)
      })

      beta_draws <- do.call(cbind, all_files)
    } else {
      beta_draws <- utils::read.table(paste0(path_folder, '/', a, '/',
                                               draws_results), head = TRUE)
    }

    variance <- utils::read.table(paste0(path_folder, '/', a, '/',
                                           model_name, '_scale.res'),
                                  head = TRUE)[1]
    dataFile <- utils::read.table(paste0(path_folder, '/', a, '/',
                                         model_name, '.data.raw'),
                                  head = TRUE)

    sigma_draws <- utils::read.table(paste0(path_folder, '/', a, '/',
                                            model_name, '_scale_sample.raw'),
                                     head = TRUE)[,2]
    spline_estimates <- NULL
    upperq_sp_estimates <- NULL
    lowerq_sp_estimates <- NULL
    if (splines){
      spline_estimates <- utils::read.table(paste0(path_folder, '/', a, '/',
                                              model_name, '_f_',
                                              name_var , '_pspline.res'),
                                            head = TRUE)[, 2:3]
      upperq_sp_estimates <-
        utils::read.table(paste0(path_folder, '/', a, '/',
                                 model_name, '_f_',
                                 name_var , '_pspline.res'),
                          head = TRUE)[, c(2, 8)]

      lowerq_sp_estimates <-
        utils::read.table(paste0(path_folder, '/', a, '/',
                                 model_name, '_f_',
                                 name_var , '_pspline.res'),
                          head = TRUE)[, c(2, 4)]
    }

    y_response <- dataFile[,'y']
    if (n_dim == 2) {
      directionX <- dataFile[, 'directionX']
    } else if (n_dim == 3){
      directionX <- dataFile[, c('directionx1', 'directionx2')]
    } else {
      directionX <- dataFile[, c('directionx1', 'directionx2', 'directionx3')]
    }

    list(fixedEffects = fixedEffects, variance = variance,
         y_response = y_response, directionX = directionX,
         fixedEffects_sd = fixedEffects_sd, beta_draws = beta_draws,
         varnames = varnames, sigma_draws = sigma_draws,
         spline_estimates = spline_estimates,
         upperq_sp_estimates = upperq_sp_estimates,
         lowerq_sp_estimates = lowerq_sp_estimates,
         lowerq = fixedEffects_lowerq,
         upperq = fixedEffects_upperq)
  })

  varnames <- results[[1]]$varnames

  directionPoint <- max(directions_ind)

  if (n_dim == 2){
    angles <- (0:(directionPoint-1))*2*pi/directionPoint
    vectorDir <- cbind(cos(angles), sin(angles))
    numbDir <- directionPoint

    orthBasis <- t(apply(vectorDir, 1, function(a){
      u_1 <- c(1,0)

      A <- cbind(a, u_1)
      qr.Q(qr(A))[, 2]
    }))
  } else if (n_dim == 3){
    x_dir <- seq(-1, 1, length = directionPoint^(1/3))
    y_dir <- seq(-1, 1, length = directionPoint^(1/3))
    z_dir <- seq(-1, 1, length = directionPoint^(1/3))

    x_y_z_grid <- expand.grid(x_dir, y_dir, z_dir)
    check_zeros <- !apply(x_y_z_grid, 1, function(a) all(a == 0))
    x_y_z_grid <- x_y_z_grid[check_zeros, ]

    vectorDir <- t(apply(x_y_z_grid, 1, function(a){
      a / sqrt(sum(a^2))
    }))
    numbDir <- nrow(vectorDir)

    orthBasis1 <- t(apply(vectorDir, 1, function(a){
      u_1 <- c(1, 0, 0)
      u_2 <- c(0, 1, 0)

      A <- cbind(a, u_1, u_2)
      qr.Q(qr(A))[, 2]
    }))

    orthBasis2 <- t(apply(vectorDir, 1, function(a){
      u_1 <- c(1, 0, 0)
      u_2 <- c(0, 1, 0)

      A <- cbind(a, u_1, u_2)
      qr.Q(qr(A))[, 3]
    }))
  } else {
    directionPoint <- directionPoint^(1/4)
    if (adaptive_dir){
      Y <- datafile_original[, response]

      x_y_z_w <- apply(Y,
                       2,
                       stats::quantile,
                       0:(directionPoint - 1)/(directionPoint - 1))

      min_xyzw <- apply(Y, 2, range)[1, ]
      max_xyzw <- apply(Y, 2, range)[2, ]

      range_xyzw <- max_xyzw - min_xyzw

      x_y_z_w_1step <- sweep(x_y_z_w, 2, min_xyzw)
      x_y_z_w_points <- (sweep(x_y_z_w_1step, 2, range_xyzw, "/") * 2) - 1

      x_y_z_w_grid <- expand.grid(x_y_z_w_points[, 1],
                                  x_y_z_w_points[, 2],
                                  x_y_z_w_points[, 3],
                                  x_y_z_w_points[, 4])

    } else {
      x_dir <- seq(-1, 1, length = directionPoint)
      y_dir <- seq(-1, 1, length = directionPoint)
      z_dir <- seq(-1, 1, length = directionPoint)
      w_dir <- seq(-1, 1, length = directionPoint)

      x_y_z_w_grid <- expand.grid(x_dir, y_dir, z_dir, w_dir)
    }

    check_zeros <- !apply(x_y_z_w_grid, 1, function(a) all(a == 0))
    x_y_z_w_grid <- x_y_z_w_grid[check_zeros, ]

    vectorDir <- t(apply(x_y_z_w_grid, 1, function(a){
      a / sqrt(sum(a^2))
    }))
    numbDir <- nrow(vectorDir)

    orthBasis1 <- t(apply(vectorDir, 1, function(a){
      u_1 <- c(1, 0, 0, 0)
      u_2 <- c(0, 1, 0, 0)
      u_3 <- c(0, 0, 1, 0)

      A <- cbind(a, u_1, u_2, u_3)
      qr.Q(qr(A))[, 2]
    }))

    orthBasis2 <- t(apply(vectorDir, 1, function(a){
      u_1 <- c(1, 0, 0, 0)
      u_2 <- c(0, 1, 0, 0)
      u_3 <- c(0, 0, 1, 0)

      A <- cbind(a, u_1, u_2, u_3)
      qr.Q(qr(A))[, 3]
    }))

    orthBasis3 <- t(apply(vectorDir, 1, function(a){
      u_1 <- c(1, 0, 0, 0)
      u_2 <- c(0, 1, 0, 0)
      u_3 <- c(0, 0, 1, 0)

      A <- cbind(a, u_1, u_2, u_3)
      qr.Q(qr(A))[, 4]
    }))
  }

  betaDifDirections_matrix <- sapply(results, function(a) a$fixedEffects)

  sdDifDirections_matrix <- sapply(results, function(a) a$fixedEffects_sd)

  lowerq_DifDirections_matrix <- sapply(results, function(a) a$lowerq)
  upperq_DifDirections_matrix <- sapply(results, function(a) a$upperq)

  posterior_sigma <- sapply(results, function(a) as.numeric(a$variance))

  sigma_draws_matrix <- sapply(results, function(a) a$sigma_draws)

  draws_matrix <- lapply(results, function(a){
    useless_columns <- which(colnames(a$beta_draws) == "intnr")
    final_draws <- a$beta_draws[, -useless_columns]
    colnames(final_draws) <- varnames
    final_draws
  })

  splines_matrix <- lapply(results, function(a) a$spline_estimates)
  upperq_splines_matrix <-lapply(results, function(a) a$upperq_sp_estimates)
  lowerq_splines_matrix <-lapply(results, function(a) a$lowerq_sp_estimates)

  organize_info <- function(object, matrix_info = TRUE){
    lapply(unique_taus, function(a){
      positions_list <- which(taus == a)
      directions_list <- directions_ind[positions_list]

      if (matrix_info){
        aux_object <- object[, positions_list]
        sapply(1:numbDir, function(b){
          aux_object[, which(directions_list == b)]
        })
      }
      else{
        aux_object <- object[positions_list]
        sapply(1:numbDir, function(b){
          aux_object[which(directions_list == b)]
        })
      }
    })
  }

  betaDifDirections <- organize_info(betaDifDirections_matrix)

  sdDifDirections <- organize_info(sdDifDirections_matrix)

  lowerq_DifDirections <- organize_info(lowerq_DifDirections_matrix)
  upperq_DifDirections <- organize_info(upperq_DifDirections_matrix)

  sigma_difDirections <- organize_info(posterior_sigma, matrix_info = FALSE)

  beta_draws_DifDirections <- organize_info(draws_matrix,
                                            matrix_info = FALSE)

  sigma_draws_DifDirections <- organize_info(sigma_draws_matrix)

  spline_estimates_DifDirections <- NULL
  upperq_spline_estimates <- NULL
  lowerq_spline_estimates <- NULL
  if (splines){
    spline_estimates_DifDirections <- organize_info(splines_matrix,
                                                    matrix_info = FALSE)
    upperq_spline_estimates <- organize_info(upperq_splines_matrix,
                                             matrix_info = FALSE)
    lowerq_spline_estimates <- organize_info(lowerq_splines_matrix,
                                             matrix_info = FALSE)
  }

  if (n_dim == 2){
    list(taus = unique_taus, betaDifDirections = betaDifDirections,
         directions = t(vectorDir), orthBases = t(orthBasis),
         sdDifDirections = sdDifDirections,
         beta_draws_DifDirections = beta_draws_DifDirections,
         sigma_difDirections = sigma_difDirections,
         sigma_draws_DifDirections = sigma_draws_DifDirections,
         spline_estimates_DifDirections = spline_estimates_DifDirections,
         lowerq_DifDirections = lowerq_DifDirections,
         upperq_DifDirections = upperq_DifDirections,
         upperq_spline_estimates = upperq_spline_estimates,
         lowerq_spline_estimates = lowerq_spline_estimates)
  } else if (n_dim == 3) {
    list(taus = unique_taus, betaDifDirections = betaDifDirections,
         directions = t(vectorDir), orthBases1 = t(orthBasis1),
         orthBases2 = t(orthBasis2),
         sdDifDirections = sdDifDirections,
         beta_draws_DifDirections = beta_draws_DifDirections,
         sigma_difDirections = sigma_difDirections,
         sigma_draws_DifDirections = sigma_draws_DifDirections,
         spline_estimates_DifDirections = spline_estimates_DifDirections,
         lowerq_DifDirections = lowerq_DifDirections,
         upperq_DifDirections = upperq_DifDirections,
         upperq_spline_estimates = upperq_spline_estimates,
         lowerq_spline_estimates = lowerq_spline_estimates)
  } else {
    list(taus = unique_taus, betaDifDirections = betaDifDirections,
         directions = t(vectorDir), orthBases1 = t(orthBasis1),
         orthBases2 = t(orthBasis2), orthBases3 = t(orthBasis3),
         sdDifDirections = sdDifDirections,
         beta_draws_DifDirections = beta_draws_DifDirections,
         sigma_difDirections = sigma_difDirections,
         sigma_draws_DifDirections = sigma_draws_DifDirections,
         spline_estimates_DifDirections = spline_estimates_DifDirections,
         lowerq_DifDirections = lowerq_DifDirections,
         upperq_DifDirections = upperq_DifDirections,
         upperq_spline_estimates = upperq_spline_estimates,
         lowerq_spline_estimates = lowerq_spline_estimates)
  }
}