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
#'  estimate. If the number of dimensions of the response variable is larger
#'  or equal to 3 then this must a matrix with all directions considered, where
#'  the rows represent each direction.
#' @param dataFile A data.frame from which to find the variables defined in the
#'  formula.
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
#' @param outfile argument to be passed to \code{bayesx.control}, in order
#'  to define a directory where all output files should be saved.
#' @param check_bayesx To check whether all calls to BayesX generated valid
#'  chain values for all models. In case there are NA values, it calls BayesX
#'  just for those models with problems. This is only considered when
#'  \code{outfile} is different than \code{NULL}.
#' @param path_bayesx When \code{check_bayes} is \code{TRUE}, the user must
#'  inform the path of BayesX in order for these new calls of the program.
#' @param ... arguments passed to \code{bayesx.control}.
#' @return A list with the chains of all parameters of interest.
#' @useDynLib baquantreg
#' @importFrom R2BayesX bayesx
#' @importFrom Formula Formula

multBayesQR <- function(response, formulaPred, directionPoint, tau = 0.5,
                        dataFile, itNum = 2000, burnin, thin = 1,
                        betaValue = NULL, sigmaValue = 1, vSampleInit = NULL,
                        priorVar = 100, hyperSigma = c(0.1, 0.1),
                        refresh = 100, bayesx = TRUE, sigmaSampling = TRUE,
                        quiet = T, tobit = FALSE, numCores = 1,
                        recordLat = FALSE, outfile = NULL,
                        check_bayesx = FALSE, path_bayesx = NULL, ...){

  n_dim <- length(response)

  if (n_dim > 3){
    numbDir <- length(directionPoint)
  } else if (length(directionPoint) > 1 & length(directionPoint) != n_dim){
    stop("Dimension of directions is different than dimension of response")
  }

  if (length(directionPoint) > 1){
    vectorDir <- directionPoint
    numbDir <- 1
  } else {
    angles <- (0:(directionPoint-1))*2*pi/directionPoint
    vectorDir <- cbind(cos(angles), sin(angles))
    numbDir <- directionPoint
  }

  objects <- list()

  objects$modelsDir <- parallel::mclapply(1:numbDir, function(a){
    if (n_dim == 2){
      if (length(directionPoint) > 1)
        u <- directionPoint
      else
        u <- vectorDir[a, ]
    } else {
      u <- directionPoint[a, ]
    }

    if (n_dim == 2){
      u_1 <- c(1, 0)

      A <- cbind(u, u_1)
      x.qr <- qr.Q(qr(A))
    }
    else {
      u_1 <- c(1, 0, 0)
      u_2 <- c(0, 1, 0)

      A <- cbind(u, u_1, u_2)
      x.qr <- qr.Q(qr(A))
    }

    Y <- dataFile[, response]

    yResp <- t(u) %*% t(Y)

    directionX <- matrix(t(x.qr[, 2:n_dim]) %*% t(Y), ncol = n_dim - 1)

    dataFile$y <- as.numeric(yResp)
    if (n_dim == 2) dataFile$directionX <- as.numeric(directionX)
    else dataFile$directionX <- directionX

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

    output$modelsTau <- lapply(tau, function(b) {
      if (bayesx){
        if (is.null(outfile)){
          result <- try(R2BayesX::bayesx(formulaUpdated, data = dataFile,
                                         iter = itNum, burnin = burnin,
                                         step = thin, method = "MCMC",
                                         family = "quantreg", quantile = b,
                                         control =
                                           R2BayesX::bayesx.control(...)))
        }
        else{
          result <- try(R2BayesX::bayesx(formulaUpdated,
                                         data = dataFile,
                                         iter = itNum, burnin = burnin,
                                         step = thin, method = "MCMC",
                                         family = "quantreg", quantile = b,
                                         control =
                                           R2BayesX::bayesx.control(...),
                                         outfile =
                                           paste0(outfile, 'dir_',
                                                  a, '_tau_', b,  '/'),
                                         dir.rm = FALSE))
        }
      }
      else {
        result <- BayesQR(tau = b, y = yResp, X = X, itNum = itNum, thin = thin,
                betaValue = betaValue, sigmaValue = sigmaValue,
                vSampleInit = vSampleInit, priorVar = priorVar,
                hyperSigma = hyperSigma, refresh = refresh,
                sigmaSampling = sigmaSampling, quiet = quiet, tobit = tobit,
                recordLat = recordLat)
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

  if (check_bayesx){
    listFolders <- list.files(outfile)

    parallel::mclapply(listFolders, function(a){
      setwd(paste0(outfile, a))

      fixed_effects_files <- grepl("_FixedEffects[0-9]+.res", list.files())
      files_results <- list.files()[fixed_effects_files]

      spline_files <- grepl("spline.res", list.files())
      splines_results <- list.files()[spline_files]
      info_splines <- utils::read.table(splines_results, head = TRUE)

      if(sum(fixed_effects_files) > 1){
        all_files <- lapply(files_results, function(aa){
          utils::read.table(aa, head = TRUE)
        })
        fixedEffects <- do.call(rbind, all_files)[, 3]
      } else {
        fixedEffects <- utils::read.table(files_results, head = TRUE)[, 3]
      }

      if(any(is.na(fixedEffect)) | any(is.na(info_splines$pmean))){
        try(system(paste(path_bayesx, 'bayesx.estim.input.prg')))
      }
    }, mc.cores = numCores)
  }


  class(objects) <- "multBQR"

  objects$method <- ifelse(bayesx, 'bayesx', 'rcpp')
  objects$response <- response
  objects$n_dim <- n_dim

  return(objects)
}