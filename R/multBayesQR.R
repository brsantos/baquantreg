#' Multiple-output Bayesian quantile regression model
#'
#' This function estimates a multiple-output Bayesian quantile regression model
#'
#' @param response Names of response variables
#' @param formulaPred a formula object, with . on the left side of a ~ operator,
#'  and the predictor terms, separated by + operators, on the right side.
#' @param tau Quantiles of interest. Default is th median, \code{tau = 0.5}.
#' @param directionPoint Either a vector with the same number of dimensions of
#'  response variable, indicating a direction, or a integer indicating the
#'  number of directions equally spaced in the unit circle one should
#'  estimate.
#' @param dataFile A data.frame from which to find the variables defined in the
#'  formula
#' @param itNum Number of iterations.
#' @param burnin Size of the initial to be discarded.
#' @param thin Thinning parameter. Default value is 1.
#' @param betaValue Initial values for the parameter beta for the continuous
#'  part.
#' @param sigmaValue Initial value for the scale parameter.
#' @param vSampleInit Initial value for the latent variables.
#' @param priorVar Value that multiplies a identity matrix in the elicition
#'  process of the prior variance of the regression parameters.
#' @param hyperSigma Vector of size containing the hyperparameters of the
#'  inverse gamma distribution for the sigma parameter of the asymmetric
#'  Laplace distribution. Default is c(0.1, 0.1), which gives a noninformative
#'  prior for sigma.
#' @param refresh Interval between printing a message during the iteration
#'  process. Default is set to 100.
#' @param bayesx If TRUE, the default, it uses bayesX software to estimate
#'  the quantile regression oarameters, which can be faster. If FALSE, it
#'  uses a Rcpp implementation of the MCMC sampler.
#' @param sigmaSampling If TRUE, the default, it will sample from the posterior
#'  distribution of the scale parameter. If FALSE, all values will be fixed to
#'  1.
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
                        priorVar = 100, hyperSigma = c(0.1, 0.1),
                        refresh = 100, bayesx = TRUE, sigmaSampling = TRUE,
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

  objects <- list()

  objects$modelsDir <- parallel::mclapply(1:numbDir, function(a){
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

    dataFile$y <- as.numeric(yResp)
    dataFile$directionX <- as.numeric(directionX)

    formulaUpdated <- stats::update(Formula::Formula(formulaPred),
                                    y ~ . + directionX)

    if(!bayesx){
      X <- stats::model.matrix(formulaUpdated, dataFile)

      if (is.null(betaValue))
        betaValue <- rep(0, dim(X)[2])
      if (is.null(vSampleInit))
        vSampleInit <- rep(1, length(yResp))
    }


    output <- list()

    output$modelsTau <- lapply(tau, function(a) {
      if (bayesx){
        # check_NA_values <- TRUE
        # while (check_NA_values){
          result <- try(R2BayesX::bayesx(formulaUpdated,
                           data = dataFile,
                           iter = itNum, burnin = burnin, step = thin,
                           method = "MCMC", family = "quantreg", quantile = a, ...))
          # check_NA_values <- FALSE

          # dimeff <- dim(result$fixed.effects)
          # if (!any(is.na(result$fixed.effects[1:dimeff[1], 1:dimeff[2]]))) check_NA_values <- FALSE
        # }
      }
      else {
        result <- BayesQR(tau = a, y = yResp, X = X, itNum = itNum, thin = thin,
                betaValue = betaValue, sigmaValue = sigmaValue, vSampleInit = vSampleInit,
                priorVar = priorVar, hyperSigma = hyperSigma,
                refresh = refresh, sigmaSampling = sigmaSampling,
                quiet = quiet,
                tobit = tobit, recordLat = recordLat, blocksV = blocksV,
                stopOrdering = stopOrdering, numOrdered = numOrdered)
      }
      result
    })

    output$tau <- tau
    output$formula <- formulaPred
    output$data <- dataFile
    output$direction <- u
    output$orthBasis = x.qr[,2]

    class(output) <- "bqr"

    return(output)
  }, mc.cores = numCores)

  class(objects) <- "multBQR"

  objects$method <- ifelse(bayesx, 'bayesx', 'rcpp')
  objects$response <- response

  return(objects)
}