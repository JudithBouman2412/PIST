## Code to recreate the concept figure manuscript

means_AUC <- determine.means()
means <- means_AUC[[2]]
AUC <- means_AUC[[1]]

# Create Concept figure

means.figure <- means[c(4,5,6)]
par(mfrow=c(1,1), mar=c(5.1,5.5,4.1,2.1))

titers <- seq(0,20,0.01)
titer.TN <- cbind(titers, dgamma(titers, shape = 1, rate = 1)/sum(dgamma(titers, shape = 1, rate = 1)))

# part A
plot(titers[1:1000], titer.TN[1:1000,2], col = "grey42", type = 'l', lwd = 5, cex.lab = 2, cex.axis = 2,
     xlab = "Quantitative test measures", ylab = 'Frequency')

cols <- c("springgreen3", "dodgerblue", "orange")

for (i in seq(1, length(means.figure),1)){
  titer.TP <- cbind(titers, dgamma(titers, shape = means.figure[i], rate = 1)/sum(dgamma(titers, shape = means.figure[i], rate = 1)))
  lines(titers[1:1000], titer.TP[1:1000,2],col = cols[i], lwd = 5)
}
legend('topright', c("Controls", "Cases -  AUC-ROC = 0.85", "Cases -  AUC-ROC = 0.9", "Cases -  AUC-ROC = 0.95")
       , col = c("grey42", cols), lwd = 5, cex = 1.5)


# part B

plot(c(0,1), c(0,1), type = 'l', lty = 2, lwd = 2.5,
     xlab = "False positive rate", ylab = "True positive rate", cex.lab=2, cex.axis =2)

for (j in seq(1, length(means.figure),1)){
  titer.TP <- cbind(titers, dgamma(titers, shape = means.figure[j], rate = 1)/sum(dgamma(titers, shape = means.figure[j], rate = 1)))

  n <- length(titer.TN[,1])

  FPR = rep(0, n)
  TPR = rep(0, n)

  for (i in seq(1,n,1)){
    TN <- sum(titer.TN[1:i,2])
    TP <- sum(titer.TP[(i:length(titer.TN[,1])),2])
    FN <- sum(titer.TP[1:i,2])
    FP <- sum(titer.TN[(i:length(titer.TN[,1])),2])
    TPR[i] <- TP/(TP+FN)
    FPR[i] <- FP/(TN+FP)
  }

  TPR <- rev(TPR)
  FPR <- rev(FPR)

  lines(FPR, TPR, col = cols[j], lwd = 5)
}

# Part C
plot(titers[1:1000], titer.TN[1:1000,2], col = "grey42", type = 'l', lwd = 5, cex.lab = 2, cex.axis = 2,
     xlab = "Raw serological test output", ylab = 'Frequency')
titer.TP <- cbind(titers, dgamma(titers, shape = means.figure[3], rate = 1)/sum( dgamma(titers, shape = means.figure[3], rate = 1)))
lines(titers[1:1000], titer.TP[1:1000,2], col = cols[3], lwd = 5)

# add threshold values:

youden <- classical.threshold(titer.TN = titer.TN, titer.TP = titer.TP)[1]
specificity <- classical.threshold(titer.TN = titer.TN, titer.TP = titer.TP, method = 'specificity')[1]

linecolors = c("turquoise3","maroon")

lines(rep(youden,2) , c(0,0.01), lty = 3, lwd = 3, col = linecolors[1])
lines(rep(specificity,2) , c(0,0.01), lty = 4, lwd = 3, col = linecolors[2])


legend('topright', c("Controls", "Cases -  AUC-ROC = 0.95", "Max Youden", "High specificity")
       , col = c("grey42", cols[3], linecolors), lwd = c(5,5,3,3,3,3) , cex = 1.5, lty=c(1,1,3,4,5,6))

# Part D

plot(c(0,1), c(0,1), type = 'l', lty = 2, lwd = 2.5,
     xlab = "False positive rate", ylab = "True positive rate", cex.lab=2, cex.axis =2)
titer.TP <- cbind(titers, dgamma(titers, shape = means.figure[3], rate = 1))

n <- length(titer.TN[,1])

FPR = rep(0, n)
TPR = rep(0, n)

for (i in seq(1,n,1)){
  TN <- sum(titer.TN[1:i,2])
  TP <- sum(titer.TP[(i:length(titer.TN[,1])),2])
  FN <- sum(titer.TP[1:i,2])
  FP <- sum(titer.TN[(i:length(titer.TN[,1])),2])
  TPR[i] <- TP/(TP+FN)
  FPR[i] <- FP/(TN+FP)
}

TPR <- rev(TPR)
FPR <- rev(FPR)

lines(FPR, TPR, col = cols[3], lwd = 5)

# add points for the various thresholds
k <- seq(1,n,1)[abs(titers-youden) == min(abs(titers-youden))]
TN <- sum(titer.TN[1:k,2])
TP <- sum(titer.TP[(k:length(titer.TN[,1])),2])
FN <- sum(titer.TP[1:k,2])
FP <- sum(titer.TN[(k:length(titer.TN[,1])),2])
TPR.youden <- TP/(TP+FN)
FPR.youden <- FP/(TN+FP)

points(FPR.youden, TPR.youden, lwd = 5, col = linecolors[1])

k <- seq(1,n,1)[abs(titers-specificity) == min(abs(titers-specificity))]
TN <- sum(titer.TN[1:k,2])
TP <- sum(titer.TP[(k:length(titer.TN[,1])),2])
FN <- sum(titer.TP[1:k,2])
FP <- sum(titer.TN[(k:length(titer.TN[,1])),2])
TPR.specificity <- TP/(TP+FN)
FPR.specificity <- FP/(TN+FP)

points(FPR.specificity, TPR.specificity, lwd = 5, col = linecolors[2])

legend('bottomright', c("Controls", "Cases -  AUC-ROC = 0.95", "Max Youden", "High specificity")
       , col = c("grey42", cols[3], linecolors), lwd = c(5,5,5,5) , cex = 1.5, lty=c(1,1,0,0), pch=c(29,29,1,1))

