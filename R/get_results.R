## Organizing the results in a folder

get_results <- function(path_folder, model_name = "bayesx.estim"){
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
  taus <- as.numeric(ifelse(length_tau == 7, substr(folders, start_positions +4,
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
      varnames <- do.call(rbind, all_files)[, 2]
    } else {
      info <- utils::read.table(paste0(path_folder, '/', a, '/',
                                       files_results), head = TRUE)
      fixedEffects <- info[, 3]
      fixedEffects_sd <- info[, 4]
      varnames <- info[, 2]
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

    sigma <- utils::read.table(paste0(path_folder, '/', a, '/',
                                         model_name, '_scale.res'),
                                  head = TRUE)[1]
    dataFile <- utils::read.table(paste0(path_folder, '/', a, '/',
                                         model_name, '.data.raw'),
                                  head = TRUE)

    y_response <- dataFile[,'y']
    directionX <- dataFile[, 'directionX']

    list(fixedEffects = fixedEffects, sigma = sigma,
         y_response = y_response, directionX = directionX,
         fixedEffects_sd = fixedEffects_sd, beta_draws = beta_draws,
         varnames = varnames)
  })

  varnames <- results[[1]]$varnames

  directionPoint <- max(directions_ind)

  angles <- (0:(directionPoint-1))*2*pi/directionPoint
  vectorDir <- cbind(cos(angles), sin(angles))
  numbDir <- directionPoint

  orthBasis <- t(apply(vectorDir, 1, function(a){
    u_1 <- c(1,0)

    A <- cbind(a, u_1)
    qr.Q(qr(A))[, 2]
  }))

  index <- sample(1:numbDir, 1)
  u <- vectorDir[index, ]
  uu <- orthBasis[index, ]

  y_resp <- results[[which(directions_ind == index)[1]]]$y_response
  xDirec <- results[[which(directions_ind == index)[1]]]$directionX

  Y <- solve(rbind(u, uu)) %*% rbind(y_resp, xDirec)

  betaDifDirections_matrix <- sapply(1:length(taus), function(a){
    results[[a]]$fixedEffects
  })

  sdDifDirections_matrix <- sapply(1:length(taus), function(a){
    results[[a]]$fixedEffects_sd
  })

  draws_matrix <- lapply(1:length(taus), function(a){
    useless_columns <- which(colnames(results[[a]]$beta_draws) == "intnr")
    final_draws <- results[[1]]$beta_draws[,-useless_columns]
    colnames(final_draws) <- varnames
    final_draws
  })

  betaDifDirections <- lapply(unique_taus, function(a){
    positions_list <- which(taus == a)
    betaSubmatrix <- betaDifDirections_matrix[, positions_list]
    directions_list <- directions_ind[positions_list]
    sapply(1:numbDir, function(b){
      betaSubmatrix[, which(directions_list == b)]
    })
  })

  sdDifDirections <- lapply(unique_taus, function(a){
    positions_list <- which(taus == a)
    sdSubmatrix <- sdDifDirections_matrix[, positions_list]
    directions_list <- directions_ind[positions_list]
    sapply(1:numbDir, function(b){
      sdSubmatrix[, which(directions_list == b)]
    })
  })

  draws_DifDirections <- lapply(unique_taus, function(a){
    positions_list <- which(taus == a)
    draws_sublist <- draws_matrix[positions_list]
    directions_list <- directions_ind[positions_list]
    lapply(1:numbDir, function(b){
      draws_sublist[which(directions_list == b)]
    })
  })

  list(Y = Y, taus = unique_taus, betaDifDirections = betaDifDirections,
       directions = t(vectorDir), orthBases = t(orthBasis), sdDifDirections =
         sdDifDirections, draws_DifDirections = draws_DifDirections)
}