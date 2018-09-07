## Calculate quantile region.

calculateQuantileRegion <- function(modelEstimates, y, numberDirections = dim(modelEstimates)[2], 
                                    ngridpoints = 100, xValue){
  
  angles <- (0:(numberDirections-1))*2*pi/numberDirections
  vectorDir <- cbind(cos(angles), sin(angles))
  
  orthBasis <- t(apply(vectorDir, 1, function(a){
    u_1 <- c(1,0)
    
    A <- cbind(a, u_1)
    x.qr <- qr.Q(qr(A))[, 2]
  }))
  
  y1range <- range(y[, 1])
  y2range <- range(y[, 2])
  
  seqY1 <- seq(y1range[1], y1range[2], length.out = ngridpoints)
  seqY2 <- seq(y2range[1], y2range[2], length.out = ngridpoints)
  
  checkPoints(seqY1, seqY2, as.numeric(tail(modelEstimates, 1)), 
              vectorDir, orthBasis, modelEstimates[-dim(modelEstimates)[1], ], 
              xValue)
}

