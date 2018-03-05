plot.summary.multBQR <- function(allEstimates){

  getResults <- function(varName, summary_multBQR){
    position <- summary_multBQR[[1]]$BetaPosterior[[1]][,1] == varName
    numDirections <- length(summary_multBQR)
    taus <- sapply(summary_multBQR[[1]], names)$BetaPosterior

    estimates <- as.numeric(sapply(summary_multBQR, function(a){
      sapply(a$BetaPosterior, function(b){
        b[position, 2:4]
      })
    }))

    data.frame(estimates = estimates,
               colorScheme = estimates > 0,
               taus = rep(taus, each = 3, times = numDirections),
               directionInd = rep(1:numDirections, each = length(taus)*3),
               typeEstimate = rep(c('coef', 'lower', 'upper'), times = numDirections*length(taus)),
               typeCoded = rep(c(1,2,2), times = numDirections*length(taus)))
  }

  namesVariables <- allEstimates[[1]]$BetaPosterior[[1]]$variable
  namesVariables <- namesVariables[namesVariables != 'directionX']

  lapply(namesVariables, function(b){
    estimatesVar <- getResults(b, allEstimates)

    numDirection <- length(allEstimates)
    angles <- (0:(numDirection-1))*2*pi/numDirection

    lapply(levels(estimatesVar$taus), function(a){
      subData <- subset(estimatesVar, taus == a)

      maxValue <- max(abs(subData$estimates))*1.05

      subData$newEstimate <- abs(subData$estimates)/maxValue

      subData$angles <- rep(angles, each = 3)

      subData$xCoord <- subData$newEstimate * cos(subData$angles)
      subData$yCoord <- subData$newEstimate * sin(subData$angles)

      cosines <- cos(c((0:(100-1))*2*pi/100, 0))
      sines <- sin(c((0:(100-1))*2*pi/100, 0))

      ggplot() + theme_bw() + ggtitle(label = paste("Variable = ", b, ", ", a, sep = "")) +
        geom_segment(data = data.frame(x = subset(subData, typeEstimate == 'lower')$xCoord,
                                       y = subset(subData, typeEstimate == 'lower')$yCoord,
                                       xend = subset(subData, typeEstimate == 'upper')$xCoord,
                                       yend = subset(subData, typeEstimate == 'upper')$yCoord),
                     aes(x = x, y = y, xend = xend, yend = yend), size = 0.5, color = 'grey75') +
        geom_point(data = subData,
                   aes(x = xCoord, y = yCoord, shape = factor(typeCoded), color = colorScheme), size = 2)  +
        theme(legend.position = 'none') + ylim(c(-1,1)) + xlim(c(-1,1)) + xlab("") + ylab("") +
        theme(axis.ticks = element_blank(), axis.text = element_blank()) +
        geom_text(data = data.frame(x = 0, y = 0, label = 0), aes(x = x, y = y, label = label)) +
        geom_text(data = data.frame(x = 1, y = 0, label = round(maxValue, 2)), aes(x = x, y = y, label = label))
    })
  })
}