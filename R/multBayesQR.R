#' Multiple-output Bayesian quantile regression model
#'
#' This function estimates a multiple-output Bayesian quantile regression model
#'
#' @param response Names of response variables
#' @param formulaPred a formula object, with . on the left side of a ~ operator,
#'  and the predictor terms, separated by + operators, on the right.
#' @param tau Quantile of interest.
#' @param directionPoint Either a vector with the same number of dimensions of
#'  response variable, indicating a direction, or a integer indicating the
#'  number of directions equally spaced in the unit circle one should
#'  estimate.
#' @param dataFile A data.frame from which to find the variables defined in the
#' formula
#' @param itNum Number of iterations.
#' @param burnin Size of the initial to be discarded.
#' @param thin Thinning parameter. Default value is 1.
#' @param betaValue Initial values for the parameter beta for the continuous
#' part.
#' @param sigmaValue Initial value for the scale parameter.
#' @param vSampleInit Initial value for the latent variables.
#' @param priorVar Value that multiplies a identity matrix in the elicition
#' process of the prior variance of the regression parameters.
#' @param refresh Interval between printing a message during the iteration
#' process. Default is set to 100.
#' @param quiet If TRUE, the default, it does not print messages to check if
#'  the MCMC is actually updating. If FALSE, it will use the value of refresh
#'  to print messages to control the iteration process.
#' @param tobit If TRUE, it will input the censored value for all observations
#'  with y = 0, according to the model. If FALSE, the default, it will estimate
#'  the parameter without this inputation process.
#' @param numCores The number of cores that could be used for estimating
#'  parallel models when more than one direction is considered.
#' @param recordLat If TRUE, it will keep the Markov chain samples for the
#'  latent variable. Default is FALSE.
#' @param blocksV Number of blocks used to sample in the posterior distribution
#'  of the latent variable. If 0, then blocking is not used and all latent
#'  observations are sampled from. Default value is 0.
#' @param stopOrdering If TRUE, it will stop ordering the weighted residuals
#'  in order to update the states of the latent variables, and will consider
#'  the ordering of some particular state of the chain; if FALSE, for every
#'  iteration of the MCMC procedure, it will keep reordering these residual
#'  terms. Default is FALSE.
#' @param numOrdered The number of iterations that will be used to order
#'  the weighted residuals needed for the update of the posterior
#'  distribution of the latent variables. Default is half the size of
#'  the MCMC chain.
#' @param ... Options to be passed to \code{bayesx} call.
#' @return A list with the chains of all parameters of interest.
#' @useDynLib baquantreg
#' @importFrom R2BayesX bayesx
#' @importFrom Formula Formula

multBayesQR <- function(response, formulaPred, directionPoint, tau = 0.5, dataFile, itNum = 2000,
                        burnin, thin = 1,
                        betaValue = NULL, sigmaValue = 1, vSampleInit = NULL,
                        priorVar = 100, refresh = 100,
                        quiet = T, tobit = FALSE, numCores = 1, recordLat = FALSE,
                        blocksV = 0, stopOrdering = FALSE, numOrdered = itNum/2, ...){

  if (length(directionPoint) > 1){
    vectorDir <- directionPoint
    numbDir <- 1
  }
  else{
    angles <- (0:(directionPoint-1))*2*pi/directionPoint
    vectorDir <- cbind(cos(angles), sin(angles))
    numbDir <- directionPoint
  }

  objects <- parallel::mclapply(1:numbDir, function(a){
    if (length(directionPoint) > 1) u <- directionPoint
    else u <- vectorDir[a,]

    u_1 <- c(1,0)

    A <- cbind(u, u_1)
    x.qr <- qr.Q(qr(A))

    # Y <- stats::model.extract(stats::model.frame(formula, dataFile),
    #                               "response")

    Y <- dataFile[, response]

    yResp <- t(u) %*% t(Y)

    directionX <- matrix(t(x.qr[,2]) %*% t(Y), ncol = 1)

    # X <- stats::model.matrix(formula, dataFile)
    #
    # X <- cbind(X, directionX)
    #
    # if (is.null(betaValue))
    #   betaValue <- rep(0, dim(X)[2])
    # if (is.null(vSampleInit))
    #   vSampleInit <- rep(1, length(yResp))

    output <- list()

    # output$chains <- lapply(tau, function(a) {
    #   BayesQR(tau = a, y = yResp, X = X, itNum = itNum, thin = thin,
    #           betaValue = betaValue, sigmaValue = sigmaValue, vSampleInit = vSampleInit,
    #           priorVar = priorVar, refresh = refresh, quiet = quiet,
    #           tobit = tobit, recordLat = recordLat, blocksV = blocksV,
    #           stopOrdering = stopOrdering, numOrdered = numOrdered)
    # })

    dataFile$y <- as.numeric(yResp)
    dataFile$directionX <- as.numeric(directionX)

    output <- lapply(tau, function(a) {
      R2BayesX::bayesx(stats::update(Formula::Formula(formulaPred),
                                               y ~ . + directionX),
                       data = dataFile,
                       iter = itNum, burnin = burnin, step = thin,
                       method = "MCMC", family = "quantreg", quantile = a, ...)

      # BayesQR(tau = a, y = yResp, X = X, itNum = itNum, thin = thin,
      #         betaValue = betaValue, sigmaValue = sigmaValue, vSampleInit = vSampleInit,
      #         priorVar = priorVar, refresh = refresh, quiet = quiet,
      #         tobit = tobit, recordLat = recordLat, blocksV = blocksV,
      #         stopOrdering = stopOrdering, numOrdered = numOrdered)
    })

    # output$tau <- tau
    # output$formula <- formula
    # output$data <- dataFile
    output$direction <- u
    output$orthBasis = x.qr[,2]
    # class(output) <- "bqr"
    return(output)
  }, mc.cores = numCores)

  class(objects) <- "multBQR"

  return(objects)
}