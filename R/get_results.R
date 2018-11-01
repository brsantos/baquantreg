## Organizing the results in a folder

get_results <- function(path_folder, model_name = "bayesx.estim"){
  folders <- list.files(path_folder)

  ## Considering at most 999 directions
  length_directions <- attr(regexpr("dir_[0-9]+_tau", folders), 'match.length')
  directions_ind <- as.numeric(ifelse(length_directions == 10, substr(folders, 5, 6),
                                      ifelse(length_directions == 9, substr(folders, 5, 5), substr(folders, 5, 7))))

  # Considering taus with at most two decimal places
  start_positions <- regexpr("tau_0.[0-9]+", folders)
  length_tau <- attr(regexpr("tau_0.[0-9]+", folders), 'match.length')
  taus <- as.numeric(ifelse(length_tau == 7, substr(folders, start_positions+4, start_positions+6),
                            substr(folders, start_positions+4, start_positions+7)))
  unique_taus <- unique(taus)

  ## Going through all folders and getting all estimates.
  results <- lapply(folders, function(a){
    fixedEffects1 <- utils::read.table(paste0(path_folder, '/', a, '/', model_name, '_FixedEffects1.res'), head = TRUE)
    fixedEffects2 <- utils::read.table(paste0(path_folder, '/', a, '/', model_name, '_FixedEffects2.res'), head = TRUE)
    fixedEffects <- rbind(fixedEffects1, fixedEffects2)[, 3]
    variance <- utils::read.table(paste0(path_folder, '/', a, '/', model_name, '_scale.res'), head = TRUE)
    dataFile <- utils::read.table(paste0(path_folder, '/', a, '/', model_name, '.data.raw'), head = TRUE)

    y_response <- dataFile[,'y']
    directionX <- dataFile[, 'directionX']

    list(fixedEffects = fixedEffects, variance = variance,
         y_response = y_response, directionX = directionX)
  })

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

  betaDifDirections <- lapply(unique_taus, function(a){
    positions_list <- which(taus == a)
    betaSubmatrix <- betaDifDirections_matrix[, positions_list]
    directions_list <- directions_ind[positions_list]
    sapply(1:numbDir, function(b){
      betaSubmatrix[, which(directions_list == b)]
    })
  })

  list(Y = Y, taus = unique_taus, betaDifDirections = betaDifDirections,
       directions = t(vectorDir), orthBases = t(orthBasis))
}