#### the functions in this file are used to evaluate titer data obtained from antibody tests #####


#' get.st
#'
#' @param fit result of the optim function
#'
#' @return matrix with the values and upper and lower bands for the
#' @export
#'
#' @examples
get.st <- function(fit){
  # function to get the confidence interval on g0 and g1 given the output of optim
  fisher_info<-solve(-fit$hessian)
  prop_sigma<-sqrt(diag(fisher_info))

  return(prop_sigma)
}

#' ll.titer
#'
#' Calculates the log likelihood of observing the quantitative test measures (titer.data) for a
#' prevalence given the distribution of quantitive test measures for the controls and cases
#'
#' @param prev the seroprevalence for which the loglikelihood is calculated
#' @param mean the shape parameter of the distribution of the case sera
#' @param scale.cases the scale parameter of the distirbution of the case sera
#' @param shape.controls the shape parameter of the distribution of the control sera
#' @param scale.controls the scale parameter of the distribution of the control sera
#' @param data.serosurvey the reuslts of the serosurvey
#'
#' @return loglikelihood
#' @export
#'
#' @examples
ll.titer <- function( prev, mean, scale.cases, shape.controls, scale.controls, data.serosurvey){

  # function calculating the likelihood of observing this data given the prevalence
  l <- prev * dgamma( data.serosurvey, shape = mean, scale = scale.cases ) +
    (1-prev) * dgamma( data.serosurvey, shape = shape.controls, scale = scale.controls ) #option log in dgamma

  ll <- sum(log(l)) # check log1p function

  # check this function by plotting the likelihood over various values of prev.

  return(ll)
}

#' analyze.serosurvey
#'
#' analyze quantitative test measures with the likelihood-based method or one of the cutoff-based methods
#'
#' @param data.serosurvey data of all observed quantitative test measures
#' @param mean the mean of the distribution of cases
#' @param type method used either 'naive' (cutoff-based) or 'likelihood'
#' @param threshold in case of 'naive' the threshold to assign someone a positive seroconversion status
#' @param mean the mean of the distribution of case sera
#' @param scale.cases the scale parameter of the distirbution of the case sera
#' @param shape.controls the shape parameter of the distribution of the control sera
#' @param scale.controls the scale parameter of the distribution of the control sera
#' @param correct TRUE if the estimate should be corrected for the sensitivity and specificity
#'
#' @return a vector with the estimate of the seroprevalence and the standard deviation in the estimate
#' @export
#'
#' @examples
analyze.serosurvey <- function(data.serosurvey = data.serosurvey, mean = 4, scale.cases = 1, shape.controls = 1, scale.controls = 1,
                          type = 'likelihood', threshold = NULL, TNR = NULL, TPR = NULL, correct = FALSE){

  if (type == 'likelihood'){
    # initial guess of prevalence
    prev_init <- sum(data.serosurvey>(max(data.serosurvey)/3))/length(data.serosurvey)

    fit <- try(optim(par = c(prev = prev_init), fn = function(x){ ll.titer( prev = x[1], mean, scale.cases, shape.controls, scale.controls, data.serosurvey) },
                     control = list(fnscale = -1), method = 'Brent', lower = 0, upper = 1, hessian = TRUE)) # try MLE2 instead of optim
    prev_est <- fit$par
    #optain CI interval for the fit
    sds <- get.st(fit)
  } else if ((type == 'naive') & correct){
    prev_est <- (sum(data.serosurvey>threshold)/length(data.serosurvey) + TNR - 1  ) / ( TPR + TNR - 1 )
    sds <- NA
  } else if ((type == 'naive') & !correct){
    prev_est <- (sum(data.serosurvey>threshold)/length(data.serosurvey) )
    sds <- NA
  }

  return(c(prev_est, sds))
}

#' sim.serosurvey
#'
#' @param N_tests number of individuals enrolled in the serosurvey
#' @param mean shape parameter of the case sera
#' @param prev the seroprevalence
#'
#' @return vector with the quantitative test values of each individual
#' @export
#'
#' @examples
sim.serosurvey <- function(N_tests = 1e4, mean = 4, prev = 0.08){
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


#' classical.threshold
#'
#' Calculates the threshold for either the maximal youden or the high specificity method
#'
#' @param mean the mean of the distribution of case sera
#' @param method 'youden' or 'specificity'
#'
#' @return a vector with the threshold value, the true positive rate and the true negative rate
#' @export
#'
#' @examples
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

#' AUC.ROC
#'
#' Calculates the area under the curve of a ROC curve.
#'
#' @param mean the shape parameter of the distribution of case sera
#'
#' @return the area under the ROC-curve
#' @export
#'
#' @examples
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


#' determine.means
#'
#' Calculate the required means for the distribution of the cases to obtain the desired AUC-ROC curves
#'
#' @param desired.AUC the AUC values that should be the result of the distribution of case sera
#'
#' @return list with as first element the true AUC-ROC values and as second element the means corresponding to those AUC values
#' @export
#'
#' @examples
determine.means <- function(desired.AUC = c(0.7, 0.75,0.8,0.85,0.9,0.95,0.975, 1) ){

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


#' create.figure
#'
#' Make figure with inferred prevalences for all 4 methods and various AUC-ROC values.
#'
#' @param N_tests number of individuals tested in each simulation
#' @param N number of simulations performed for each AUC value
#' @param means the means of the titer distribution that match the AUC values
#' @param AUC the AUC-ROC values that should be tested
#' @param true_prev the true prevalence
#'
#' @return the figure
#' @export
#'
#' @examples
create.figure <- function(N_tests = 1000, N = 50, means, AUC, true_prev=0.08){

  #create matrixes for results
  results_value <- matrix(0, nrow = length(means), ncol = N)
  results_sd <- matrix(0, nrow = length(means), ncol = N)

  results_naive_classical <- matrix(0, nrow = length(means), ncol = N)
  results_naive_spec <- matrix(0, nrow = length(means), ncol = N)

  results_naive_classical_uncor <- matrix(0, nrow = length(means), ncol = N)
  results_naive_spec_uncor <- matrix(0, nrow = length(means), ncol = N)

  #simulate results for likelihood method
  for (i in seq(1, length(means), 1)){

    for (j in seq(1,N,1)){
      # simulate data given those distributions
      data.serosurvey <- sim.serosurvey( N_tests = N_tests,  mean = means[i], prev=true_prev)

      # analyze data and estimate prev
      output <- analyze.serosurvey(data.serosurvey = data.serosurvey, mean = means[i], type = 'likelihood')

      # what is the probability that our estimated prevalence is equal to the true prevalence?

      results_value[i,j] <- output[[1]]
      results_sd[i,j] <- output[[2]]

      # use as a threshold the crossing between both curves --> different choises would have a different effect!
      threshold.clas = classical.threshold(mean = means[i])
      threshold.spec = classical.threshold(mean = means[i], method = 'specificity')

      # analyze data and estimate prev
      output.clas <- analyze.serosurvey(data.serosurvey = data.serosurvey, mean = means[i], type = 'naive',
                             threshold = threshold.clas[1], TNR = threshold.clas[3], TPR = threshold.clas[2], correct = TRUE)
      output.spec <- analyze.serosurvey(data.serosurvey = data.serosurvey, mean = means[i], type = 'naive',
                             threshold = threshold.spec[1], TNR = threshold.spec[3], TPR = threshold.spec[2], correct = TRUE)

      # what is the probability that our estimated prevalence is equal to the true prevalence?
      results_naive_classical[i,j] <- output.clas[[1]]
      results_naive_spec[i,j] <- output.spec[[1]]

      output.clas_uncor <- analyze.serosurvey(data.serosurvey = data.serosurvey, mean = means[i], type = 'naive',
                                   threshold = threshold.clas[1], TNR = threshold.clas[3], TPR = threshold.clas[2], correct = FALSE)
      output.spec_uncor <- analyze.serosurvey(data.serosurvey = data.serosurvey, mean = means[i], type = 'naive',
                                   threshold = threshold.spec[1], TNR = threshold.spec[3], TPR = threshold.spec[2], correct = FALSE)

      # what is the probability that our estimated prevalence is equal to the true prevalence?
      results_naive_classical_uncor[i,j] <- output.clas_uncor[[1]]
      results_naive_spec_uncor[i,j] <- output.spec_uncor[[1]]
    }

  }

  # Create dataframe
  results.ll <- matrix(0, nrow = length(means)*N, ncol = 2)
  for (i in seq(1, length(means), 1)){
    for (j in seq(1, N,1 )){
      results.ll[ ((i-1)*N + j) ,1] <- results_value[i,j]
      results.ll[ ((i-1)*N + j) ,2] <- round(AUC[i], 2)
    }
  }

  results.ll <- as.data.frame(results.ll)
  colnames( results.ll ) <- c('Value','AUC')
  results.ll[,2] <- as.factor(results.ll[,2])
  results.ll[,3] <- "Likelihood"

  # Create dataframes
  results.ll.clas <- matrix(0, nrow = length(means)*N, ncol = 2)
  for (i in seq(1, length(means), 1)){
    for (j in seq(1, N,1 )){
      results.ll.clas[ ((i-1)*N + j) ,1] <- results_naive_classical[i,j]
      results.ll.clas[ ((i-1)*N + j) ,2] <- round(AUC[i], 2)
    }
  }

  results.ll.clas <- as.data.frame(results.ll.clas)
  colnames( results.ll.clas ) <- c('Value','AUC')
  results.ll.clas[,2] <- as.factor(results.ll.clas[,2])
  results.ll.clas[,3] <- "Max Youden Corrected"

  results.ll.spec <- matrix(0, nrow = length(means)*N, ncol = 2)
  for (i in seq(1, length(means), 1)){
    for (j in seq(1, N,1 )){
      results.ll.spec[ ((i-1)*N + j) ,1] <- results_naive_spec[i,j]
      results.ll.spec[ ((i-1)*N + j) ,2] <- round(AUC[i], 2)
    }
  }

  results.ll.spec <- as.data.frame(results.ll.spec)
  colnames( results.ll.spec ) <- c('Value','AUC')
  results.ll.spec[,2] <- as.factor(results.ll.spec[,2])
  results.ll.spec[,3] <- "High specificity Corrected"

  # Create dataframes
  results.ll.clas_uncor <- matrix(0, nrow = length(means)*N, ncol = 2)
  for (i in seq(1, length(means), 1)){
    for (j in seq(1, N,1 )){
      results.ll.clas_uncor[ ((i-1)*N + j) ,1] <- results_naive_classical_uncor[i,j]
      results.ll.clas_uncor[ ((i-1)*N + j) ,2] <- round(AUC[i], 2)
    }
  }

  results.ll.clas_uncor <- as.data.frame(results.ll.clas_uncor)
  colnames( results.ll.clas_uncor ) <- c('Value','AUC')
  results.ll.clas_uncor[,2] <- as.factor(results.ll.clas_uncor[,2])
  results.ll.clas_uncor[,3] <- "Max Youden"

  results.ll.spec_uncor <- matrix(0, nrow = length(means)*N, ncol = 2)
  for (i in seq(1, length(means), 1)){
    for (j in seq(1, N,1 )){
      results.ll.spec_uncor[ ((i-1)*N + j) ,1] <- results_naive_spec_uncor[i,j]
      results.ll.spec_uncor[ ((i-1)*N + j) ,2] <- round(AUC[i], 2)
    }
  }

  results.ll.spec_uncor <- as.data.frame(results.ll.spec_uncor)
  colnames( results.ll.spec_uncor ) <- c('Value','AUC')
  results.ll.spec_uncor[,2] <- as.factor(results.ll.spec_uncor[,2])
  results.ll.spec_uncor[,3] <- "High specificity"


  ll.all <- rbind(results.ll, results.ll.clas, results.ll.spec, results.ll.clas_uncor, results.ll.spec_uncor)
  names(ll.all)[3] <- 'Type'
  ll.all[,3] <- factor(ll.all[,3], levels = c(  "Max Youden", "Max Youden Corrected","High specificity","High specificity Corrected",  "Likelihood" ))

  plot.out <- ggplot() +
    geom_violin(data=ll.all, aes(x = AUC, y = Value, fill=Type,colour=Type), width = 1, alpha = .2, size=0.3)+
    theme(axis.text = element_text(size=26),axis.title = element_text(size=26), text = element_text(size=26),legend.position="bottom",
          panel.grid.major =element_line(colour = "grey"), panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(), legend.background = element_blank(),plot.title = element_text(hjust = 0.5))+
    geom_hline(aes(yintercept = true_prev, color = "True prevalence"), size = 1 )+
    ylab("Estimated seroprevalence") + xlab("Area under ROC-curve (AUC - ROC)") +
    scale_colour_manual("", values = c("High specificity"="violetred3","Max Youden" = "dodgerblue" ,
                                       "High specificity Corrected"="orange", "Max Youden Corrected" = "cyan", "Likelihood"="springgreen3", "True prevalence" = "red" ),
                        aesthetics = c("colour","fill"))
  return(plot.out)

}

#' create.figure
#'
#' Make figure with inferred prevalences for all 4 methods and various AUC-ROC values.
#'
#' @param N_tests number of individuals tested in each simulation
#' @param N number of simulations performed for each AUC value
#' @param means the means of the titer distribution that match the AUC values
#' @param AUC the AUC-ROC values that should be tested
#' @param true_prev the true prevalence
#'
#' @return the figure
#' @export
#'
#' @examples
create.figure.bootstrap <- function(N_tests = 10000, N = 200, means, AUC, true_prev = 0.08){

  #create matrixes for results
  results_value.bootstrap <- matrix(0, nrow = length(means), ncol = N)
  results_sd.bootstrap <- matrix(0, nrow = length(means), ncol = N)

  results_naive_classical.bootstrap <- matrix(0, nrow = length(means), ncol = N)
  results_naive_spec.bootstrap <- matrix(0, nrow = length(means), ncol = N)

  #simulate results for likelihood method
  for (i in seq(1, length(means), 1)){
    # simulate data given those distributions
    data.serosurvey.original <- sim.serosurvey( N_tests = N_tests,  mean = means[i], prev=true_prev)
    # use as a threshold the crossing between both curves --> different choises would have a different effect!
    threshold.clas = classical.threshold(mean = means[i])
    threshold.spec = classical.threshold(mean = means[i], method = 'specificity')

    for (j in seq(1,N,1)){
      data.serosurvey.bootstrap <- sample(x = data.serosurvey.original, size = N_tests, replace = TRUE)

      # analyze data and estimate prev
      output <- analyze.serosurvey(data.serosurvey = data.serosurvey.bootstrap, mean = means[i], type = 'likelihood')

      # what is the probability that our estimated prevalence is equal to the true prevalence?

      results_value.bootstrap[i,j] <- output[[1]]
      results_sd.bootstrap[i,j] <- output[[2]]

      # analyze data and estimate prev
      output.clas <- analyze.serosurvey(data.serosurvey = data.serosurvey.bootstrap, mean = means[i], type = 'naive',
                                        threshold = threshold.clas[1], TNR = threshold.clas[3], TPR = threshold.clas[2], correct = TRUE)
      output.spec <- analyze.serosurvey(data.serosurvey = data.serosurvey.bootstrap, mean = means[i], type = 'naive',
                                        threshold = threshold.spec[1], TNR = threshold.spec[3], TPR = threshold.spec[2], correct = TRUE)

      # what is the probability that our estimated prevalence is equal to the true prevalence?
      results_naive_classical.bootstrap[i,j] <- output.clas[[1]]
      results_naive_spec.bootstrap[i,j] <- output.spec[[1]]

    }

  }

  # Compare with resampling the complete serosurvey every time.

  #create matrixes for results
  results_value <- matrix(0, nrow = length(means), ncol = N)
  results_sd <- matrix(0, nrow = length(means), ncol = N)

  results_naive_classical <- matrix(0, nrow = length(means), ncol = N)
  results_naive_spec <- matrix(0, nrow = length(means), ncol = N)

  #simulate results for likelihood method
  for (i in seq(1, length(means), 1)){
    # use as a threshold the crossing between both curves --> different choises would have a different effect!
    threshold.clas = classical.threshold(mean = means[i])
    threshold.spec = classical.threshold(mean = means[i], method = 'specificity')

    for (j in seq(1,N,1)){
      # simulate data given those distributions
      data.serosurvey <- sim.serosurvey( N_tests = N_tests,  mean = means[i], prev=true_prev)

      # analyze data and estimate prev
      output <- analyze.serosurvey(data.serosurvey = data.serosurvey, mean = means[i], type = 'likelihood')

      # what is the probability that our estimated prevalence is equal to the true prevalence?

      results_value[i,j] <- output[[1]]
      results_sd[i,j] <- output[[2]]

      # analyze data and estimate prev
      output.clas <- analyze.serosurvey(data.serosurvey = data.serosurvey, mean = means[i], type = 'naive',
                                        threshold = threshold.clas[1], TNR = threshold.clas[3], TPR = threshold.clas[2], correct = TRUE)
      output.spec <- analyze.serosurvey(data.serosurvey = data.serosurvey, mean = means[i], type = 'naive',
                                        threshold = threshold.spec[1], TNR = threshold.spec[3], TPR = threshold.spec[2], correct = TRUE)

      # what is the probability that our estimated prevalence is equal to the true prevalence?
      results_naive_classical[i,j] <- output.clas[[1]]
      results_naive_spec[i,j] <- output.spec[[1]]

    }

  }

  # Create dataframe
  results.ll <- matrix(0, nrow = length(means)*N, ncol = 2)
  for (i in seq(1, length(means), 1)){
    for (j in seq(1, N,1 )){
      results.ll[ ((i-1)*N + j) ,1] <- results_value[i,j] - mean(results_value[i,])
      results.ll[ ((i-1)*N + j) ,2] <- round(AUC[i], 2)
    }
  }

  results.ll <- as.data.frame(results.ll)
  colnames( results.ll ) <- c('Value','AUC')
  results.ll[,2] <- as.factor(results.ll[,2])
  results.ll[,3] <- "Likelihood"

  # Create dataframes
  results.ll.clas <- matrix(0, nrow = length(means)*N, ncol = 2)
  for (i in seq(1, length(means), 1)){
    for (j in seq(1, N,1 )){
      results.ll.clas[ ((i-1)*N + j) ,1] <- results_naive_classical[i,j] - mean(results_naive_classical[i,])
      results.ll.clas[ ((i-1)*N + j) ,2] <- round(AUC[i], 2)
    }
  }

  results.ll.clas <- as.data.frame(results.ll.clas)
  colnames( results.ll.clas ) <- c('Value','AUC')
  results.ll.clas[,2] <- as.factor(results.ll.clas[,2])
  results.ll.clas[,3] <- "Max Youden"

  results.ll.spec <- matrix(0, nrow = length(means)*N, ncol = 2)
  for (i in seq(1, length(means), 1)){
    for (j in seq(1, N,1 )){
      results.ll.spec[ ((i-1)*N + j) ,1] <- results_naive_spec[i,j] - mean(results_naive_spec[i,])
      results.ll.spec[ ((i-1)*N + j) ,2] <- round(AUC[i], 2)
    }
  }

  results.ll.spec <- as.data.frame(results.ll.spec)
  colnames( results.ll.spec ) <- c('Value','AUC')
  results.ll.spec[,2] <- as.factor(results.ll.spec[,2])
  results.ll.spec[,3] <- "High specificity"

  # Create dataframe
  results.ll.bootstrap <- matrix(0, nrow = length(means)*N, ncol = 2)
  for (i in seq(1, length(means), 1)){
    for (j in seq(1, N,1 )){
      results.ll.bootstrap[ ((i-1)*N + j) ,1] <- results_value.bootstrap[i,j] - mean(results_value.bootstrap[i,])
      results.ll.bootstrap[ ((i-1)*N + j) ,2] <- round(AUC[i], 2)
    }
  }

  results.ll.bootstrap <- as.data.frame(results.ll.bootstrap)
  colnames( results.ll.bootstrap ) <- c('Value','AUC')
  results.ll.bootstrap[,2] <- as.factor(results.ll.bootstrap[,2])
  results.ll.bootstrap[,3] <- "Likelihood bootstrap"

  # Create dataframes
  results.ll.clas.bootstrap <- matrix(0, nrow = length(means)*N, ncol = 2)
  for (i in seq(1, length(means), 1)){
    for (j in seq(1, N,1 )){
      results.ll.clas.bootstrap[ ((i-1)*N + j) ,1] <- results_naive_classical.bootstrap[i,j] - mean(results_naive_classical.bootstrap[i,])
      results.ll.clas.bootstrap[ ((i-1)*N + j) ,2] <- round(AUC[i], 2)
    }
  }

  results.ll.clas.bootstrap <- as.data.frame(results.ll.clas.bootstrap)
  colnames( results.ll.clas.bootstrap ) <- c('Value','AUC')
  results.ll.clas.bootstrap[,2] <- as.factor(results.ll.clas.bootstrap[,2])
  results.ll.clas.bootstrap[,3] <- "Max Youden bootstrap"

  results.ll.spec.bootstrap <- matrix(0, nrow = length(means)*N, ncol = 2)
  for (i in seq(1, length(means), 1)){
    for (j in seq(1, N,1 )){
      results.ll.spec.bootstrap[ ((i-1)*N + j) ,1] <- results_naive_spec.bootstrap[i,j] - mean(results_naive_spec.bootstrap[i,])
      results.ll.spec.bootstrap[ ((i-1)*N + j) ,2] <- round(AUC[i], 2)
    }
  }

  results.ll.spec.bootstrap <- as.data.frame(results.ll.spec.bootstrap)
  colnames( results.ll.spec.bootstrap ) <- c('Value','AUC')
  results.ll.spec.bootstrap[,2] <- as.factor(results.ll.spec.bootstrap[,2])
  results.ll.spec.bootstrap[,3] <- "High specificity bootstrap"

  ll.all.bootstrap <- rbind(results.ll, results.ll.clas, results.ll.spec,results.ll.bootstrap, results.ll.clas.bootstrap, results.ll.spec.bootstrap)
  names(ll.all.bootstrap)[3] <- 'Type'
  ll.all.bootstrap[,3] <- factor(ll.all.bootstrap[,3], levels = c("High specificity bootstrap","High specificity","Max Youden bootstrap","Max Youden",
                                       "Likelihood bootstrap", "Likelihood" ))
  ll.all.bootstrap[,2] <- factor(ll.all.bootstrap[,2])

  plot.out <- ggplot() +
    geom_violin(data=ll.all.bootstrap, aes(x = AUC, y = Value, fill=Type,colour=Type), width = 1, alpha = .2, size=0.3)+
    theme(axis.text = element_text(size=26),axis.title = element_text(size=26), text = element_text(size=26),legend.position="bottom",
          panel.grid.major =element_line(colour = "grey"), panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(), legend.background = element_blank(),plot.title = element_text(hjust = 0.5))+
    ylab("Deviation from mean seroprevalence est.") + xlab("Area under ROC-curve (AUC - ROC)") +
    scale_colour_manual("", values = c("High specificity bootstrap"="violetred3","High specificity"="orange","Max Youden bootstrap" = "dodgerblue" ,
                                       "Max Youden" = "darkviolet","Likelihood bootstrap"="springgreen3" ,"Likelihood"="cyan"),
                        aesthetics = c("colour","fill"))
  plot.out

  return(plot.out)

}

#' estimate.CC.distr
#'
#' Simulates the sampling of controls and cases to estimate their distribution
#'
#' @param mean mean of the distribution of case sera
#' @param sample.sizes the sample size that should be tested
#' @param name the AUC value of the tested sample sizes
#' @param N_tests
#' @param N the number of replicates
#' @param real.prev the actual prevalence
#'
#' @return a dataframe with the results of all N in silico serosurveys
#' @export
#'
#' @examples
estimate.CC.distr <- function(true.mean, sample.sizes, name="0.9",N_tests = 10000, N = 200, real.prev=0.08){
  resampling <- matrix(0, ncol=N, nrow = length(sample.sizes) )
  for (i in seq(1, length(sample.sizes),1)){
    for (j in seq(1, N, 1)){
      titer.TN.sampled <- rgamma(sample.sizes[i], shape = 1, scale = 1)
      titer.TP.sampled <- rgamma(sample.sizes[i], shape = true.mean, scale = 1)

      fit.TP <- fitdist(titer.TP.sampled, distr = "gamma")
      fit.TN <- fitdist(titer.TN.sampled, distr = "gamma", fix.arg = list(shape = 1)) # put this into methods

      # make a distribution from the sampled controls and cases
      data.serosurvey <- sim.serosurvey( N_tests = N_tests, mean = true.mean, prev=real.prev)

      # analyze data and estimate prev
      output <- try(analyze.serosurvey(data.serosurvey = data.serosurvey, mean = fit.TP[[1]][1], scale.cases = fit.TP[[1]][2],
                                  shape.controls = 1, scale.controls = fit.TN[[1]][1], type = 'likelihood'))

      # what is the probability that our estimated prevalence is equal to the true prevalence?
      resampling[i,j] <- output[[1]]
    }
  }

  # Create dataframe
  results.ll <- matrix(0, nrow = length(sample.sizes)*N, ncol = 3)
  results.ll[,3] <- name
  for (i in seq(1, length(sample.sizes), 1)){
    for (j in seq(1, N,1 )){
      results.ll[ ((i-1)*N + j) ,1] <- as.numeric(resampling[i,j])
      results.ll[ ((i-1)*N + j) ,2] <- sample.sizes[i]
    }
  }

  return(results.ll)
}
