##

multBayesQR <- function(directionPoint, tau = 0.5, Y, X, itNum = 2000, ...){

  if (length(directionPoint) > 1){
    vectorDir <- directionPoint
    numbDir <- 1
  }
  else{
    angles <- (0:(directionPoint-1))*2*pi/directionPoint
    vectorDir <- cbind(cos(angles), sin(angles))
    numbDir <- directionPoint
  }

  objects <- lapply(1:numbDir, function(a){
    if (length(directionPoint) > 1) u <- directionPoint
    else u <- vectorDir[a,]

    u_1 <- c(1,0)

    A <- cbind(u, u_1)
    x.qr <- qr.Q(qr(A))

    yResp <- t(u) %*% Y

    xDirec <- t(x.qr[,2]) %*% Y

    dataModel <- data.frame(y = as.numeric(yResp),
                            x = as.numeric(xDirec))

    model <- bayesQR(y ~ x, tau = tau, itNum = itNum, data = dataModel)
    list(direction = u, orthBasis = x.qr[,2], model = model)
  })

  list(tau = tau, objects = objects)
}