# Script to get results for Figure 2 on the cluster

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sample.size = as.numeric(args[1]) # filename to save the results
mean.i = as.numeric(args[2])
true.prev = as.numeric(args[3])

setwd('/cluster/home/jbouman/minsample')

get.CI <- function(fit){
  # function to get the confidence interval on g0 and g1 given the output of optim
  fisher_info<-solve(-fit$hessian)
  prop_sigma<-sqrt(diag(fisher_info))

  return(prop_sigma)
}

simulate.titer <- function( N_tests = 1e4, tot_pop = 8e6, mean, prev ){

  # determine how many of the tests are positive
  N_pos <- rbinom(N_tests,1,prev)

  titer.data <- rep(0, N_tests)

  # create matrix with titer values of all positive and negative values
  for (j in seq(1,N_tests,1)){
    if (N_pos[j] == 1){
      titer.data[j] <- rgamma(1, shape = mean, scale = 1)
    } else if (N_pos[j] == 0){
      titer.data[j] <- rgamma(1, shape = 1, scale = 1)
    }
  }

  return(titer.data)
}

ll.titer <- function(parms = c(prev = x[1]), extra){

  prev <- parms[1]

  mean <- extra[[1]]
  rate.TP <- extra[[2]]
  shape.TN <- extra[[3]]
  rate.TN <- extra[[4]]


  # function calculating the likelihood of observing this data given the prevalence
  l <- prev * dgamma( data.titer, shape = mean, scale = rate.TP ) +
    (1-prev) * dgamma( data.titer, shape = shape.TN, scale = rate.TN )

  ll <- sum(log(l))

  return(ll)
}

analyze.titer <- function(data.titer = data.titer, mean, rate.TP = 1, shape.TN = 1, rate.TN = 1,
                          type = 'likelihood', threshold = NULL, TNR = NULL, TPR = NULL, correct = FALSE){

  if (type == 'likelihood'){
    # initial guess of prevalence
    prev_init <- sum(data.titer>(max(data.titer)/3))/length(data.titer)

    fit <- try(optim(par = c(prev = prev_init), fn = ll.titer, extra = list(mean, rate.TP, shape.TN, rate.TN),
                     control = list(fnscale = -1), method = 'Brent', lower = 0, upper = 1, hessian = TRUE))
    prev_est <- fit$par
    #optain CI interval for the fit
    sds <- get.CI(fit)
  } else if ((type == 'naive') & correct){
    prev_est <- (sum(data.titer>threshold)/length(data.titer) + TNR - 1  ) / ( TPR + TNR - 1 )
    sds <- NA
  } else if ((type == 'naive') & !correct){
    prev_est <- (sum(data.titer>threshold)/length(data.titer) )
    sds <- NA
  }

  return(c(prev_est, sds))
}


classical.threshold <- function(mean, method = 'youden'){

  titers <- seq(0,20, 0.01)
  titer.TN <- cbind(titers, dgamma(titers, shape = 1, scale = 1))
  titer.TP <- cbind(titers, dgamma(titers, shape = mean, scale = 1))

  threshold.quality <- rep(0,length(titer.TN[,1]))
  TPR <- rep(0,length(titers))
  TNR <- rep(0,length(titers))

  for (i in seq(1, length(titers), 1)){
    TN <- sum(titer.TN[1:i,2])
    TP <- sum(titer.TP[(i:length(titer.TN[,1])),2])
    FN <- sum(titer.TP[1:i,2])
    FP <- sum(titer.TN[(i:length(titer.TN[,1])),2])
    TPR[i] <- TP/(TP+FN)
    TNR[i] <- TN/(TN+FP)
    threshold.quality[i] <- TPR[i]  +  TNR[i] - 1
  }

  if (method == 'youden'){
    threshold <- titer.TP[  threshold.quality==max(threshold.quality) ,1]
    TPR_opt <- TPR[ threshold.quality==max(threshold.quality)]
    TNR_opt <- TNR[ threshold.quality==max(threshold.quality)]
  } else if (method == 'specificity'){
    threshold <- titer.TP[    min(abs(TNR-0.99))== abs(TNR-0.99),1  ]
    TPR_opt <- TPR[ min(abs(TNR-0.99))== abs(TNR-0.99) ]
    TNR_opt <- TNR[ min(abs(TNR-0.99))== abs(TNR-0.99) ]
  }

  return(c(threshold,TPR_opt, TNR_opt))

}

AUC.ROC <- function(mean){

  titers <- seq(0,20, 0.01)
  titer.TN <- cbind(titers, dgamma(titers, shape = 1, scale = 1))
  titer.TP <- cbind(titers, dgamma(titers, shape = mean, scale = 1))

  n <- length(titers)

  FPR = rep(0, n+1)
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

  FPR[n+1] <- 1

  AUC <- sum(TPR*(FPR[2:(n+1)] - FPR[1:n]))

  return(AUC)
}


determine.means <- function(desired.AUC = c(0.7, 0.75,0.8,0.85,0.9,0.95, 0.975, 1) ){

  # distribution of true negatives
  titers <- seq(0,20,0.01)
  titer.TN <- cbind(titers, dgamma(titers, shape = 1, rate = 1))

  means <- rev(seq(1,12,length.out=500))
  AUC <- rep(0, length(means))

  for (i in seq(1, length(means),1)){
    #distribtion of true positives
    titer.TP <- cbind(titers, dgamma(titers, shape = means[i], rate = 1))

    # calculate the AUC-ROC for the TN and TP distributions
    AUC[i] <- AUC.ROC(mean = means[i])
  }

  means.final <- rep(0, length(desired.AUC))
  AUC.final <- rep(0, length(desired.AUC))

  for (i in seq(1, length(desired.AUC),1)){
    means.final[i] <- means[abs(AUC-desired.AUC[i])==min(abs(AUC-desired.AUC[i]))]
    AUC.final[i] <- AUC[abs(AUC-desired.AUC[i])==min(abs(AUC-desired.AUC[i]))]
  }

  AUC <- AUC.final
  means <- means.final

  return(list(AUC, means))
}

##### Example for using titer data #######

# example titer distributions
N = 50 # number of tests per overlap value

means <- rev(seq(1,12,length.out=1000))
AUC <- rep(0, length(means))
titers <- seq(0,20,0.01)

for (i in seq(1, length(means),1)){
  titer.TP <- cbind(titers, dgamma(titers, shape = means[i], rate = 1))

  # calculate the overlap between the TN and TP titers
  AUC[i] <- AUC.ROC(mean = means[i])
}

desired.AUC <- c(0.7,0.75,0.8,0.85,0.9,0.95, 0.975, 1)
means.final <- rep(0, length(desired.AUC))
AUC.final <- rep(0, length(desired.AUC))

for (i in seq(1, length(desired.AUC),1)){
  means.final[i] <- means[abs(AUC-desired.AUC[i])==min(abs(AUC-desired.AUC[i]))]
  AUC.final[i] <- AUC[abs(AUC-desired.AUC[i])==min(abs(AUC-desired.AUC[i]))]
}

AUC <- AUC.final[mean.i]
means <- means.final[mean.i]

# simulate data given those distributions
data.titer <- simulate.titer(N_tests = sample.size, tot_pop=8e6, mean = means, prev=true.prev)

# analyze data and estimate prev
output <- analyze.titer(data.titer = data.titer, mean = means, type = 'likelihood')

result <- t(c(AUC, sample.size, true.prev, output))

write.table(result, file = 'power_analysis_fullData_4.txt', append = TRUE,col.names = FALSE )



