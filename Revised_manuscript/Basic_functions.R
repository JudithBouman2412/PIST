
#' observed.threshold
#'
#' @param control.data data simulated from control distribution
#' @param case.data data simulated from case distribution
#' @param method method used; youden or specificity
#'
#' @return
#' @export
#'
#' @examples
observed.threshold <- function(control.data, case.data, method = 'youden'){

  titers <- seq(min(c(control.data, case.data)), max(c(control.data, case.data)),0.01)

  threshold.quality <- rep(0,length(titers))
  TPR <- rep(0,length(titers))
  TNR <- rep(0,length(titers))

  for (i in seq(1, length(titers), 1)){
    TN <- sum(control.data<titers[i])
    TP <- sum(case.data>titers[i])
    FN <- sum(case.data<titers[i])
    FP <- sum(control.data>titers[i])
    TPR[i] <- TP/(TP+FN)
    TNR[i] <- TN/(TN+FP)
    threshold.quality[i] <- TPR[i]  +  TNR[i] - 1
  }

  if (method == 'youden'){
    threshold <- titers[  threshold.quality==max(threshold.quality) ][1]
    TPR_opt <- TPR[ threshold.quality==max(threshold.quality)][1]
    TNR_opt <- TNR[ threshold.quality==max(threshold.quality)][1]
  } else if (method == 'specificity'){
    threshold <- titers[    min(abs(TNR-0.99))== abs(TNR-0.99) ][1]
    TPR_opt <- TPR[ min(abs(TNR-0.99))== abs(TNR-0.99) ][1]
    TNR_opt <- TNR[ min(abs(TNR-0.99))== abs(TNR-0.99) ][1]
  }

  return(c(threshold,TPR_opt, TNR_opt))
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
classical.threshold <- function(mean = 4, method = 'youden'){

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

#' sim.validation.data
#'
#' @param control.data data simulated from control distribution
#' @param case.data data simulated from case distribution
#'
#' @return list of case distribution and control distribution
#' @export
#'
#' @examples
sim.validation.data <- function(control.data, case.data){
  # estimates density function from control and case data

  # fit empirical distribution through it
  case.fun <- approxfun( density(case.data))
  control.fun <- approxfun( density(control.data))

  pdf.case.smooth <- function(x){
    pdf <- case.fun(x)
    pdf[is.na(pdf)]<-min(pdf, na.rm = T)
    return(pdf)
  }

  pdf.control.smooth <- function(x){
    pdf <- control.fun(x)
    pdf[is.na(pdf)]<-min(pdf, na.rm = T)
    return(pdf)
  }

  return(list(pdf.case.smooth, pdf.control.smooth))

}

#' analyze.serosurvey.emperical
#'
#' @param data.serosurvey simulated serosurvey data
#' @param pdf.control.smooth estimated distribution of control sera
#' @param pdf.case.smooth estimated distribution of case sera
#' @param type method; likelihood, youden or specificity
#' @param threshold cutoff-value for cutoff-based measures
#' @param TNR true negative rate
#' @param TPR true positive rate
#' @param correct whether the ad-hoc Rogan-Gladen should be applied
#'
#' @return
#' @export
#'
#' @examples
analyze.serosurvey.emperical <- function(data.serosurvey = data.serosurvey, pdf.control.smooth, pdf.case.smooth, type = 'likelihood',
                                         threshold = NULL, TNR = NULL, TPR = NULL, correct = FALSE){

  if (type == 'likelihood'){
    # initial guess of prevalence
    prev_init <- sum(data.serosurvey>(max(data.serosurvey)/3))/length(data.serosurvey)

    fit <- try(optim(par = c(prev = prev_init), fn = function(x){ ll.titer.emperical( prev = x[1], pdf.control.smooth, pdf.case.smooth, data.serosurvey) },
                     control = list(fnscale = -1), method = 'Brent', lower = 0, upper = 1, hessian = TRUE)) # try MLE2 instead of optim
    prev_est <- fit$par
    #optain CI interval for the fit
    return(c(prev_est, ll=fit$value))
  } else if ((type == 'youden' || type == 'specificity') & correct){
    prev_est <- (sum(data.serosurvey>threshold)/length(data.serosurvey) + TNR - 1  ) / ( TPR + TNR - 1 )
    return(c(prev_est))
  } else if ((type == 'youden' || type == 'specificity') & !correct){
    prev_est <- (sum(data.serosurvey>threshold)/length(data.serosurvey) )
    return(c(prev_est))
  }

}

#' ll.titer.emperical
#'
#' @param prev true cumulative incidence
#' @param data.serosurvey simulated serosurvey data
#' @param pdf.control.smooth estimated distribution of control sera
#' @param pdf.case.smooth estimated distribution of case sera
#'
#' @return
#' @export
#'
#' @examples
ll.titer.emperical <- function( prev, pdf.control.smooth, pdf.case.smooth , data.serosurvey){
  # likelihood function for mixture model with k = 2

  # function calculating the likelihood of observing this data given the prevalence
  l <- prev * pdf.case.smooth(data.serosurvey) +
    (1-prev) * pdf.control.smooth(data.serosurvey) #option log in dgamma

  ll <- sum(log(l)) # check log1p function

  return(ll)
}


#' bootrstrap.analyze
#'
#' Create bootstrap estimate and CI interval for all non-bayesian analyses
#'
#' @param data.serosurvey simulated serosurvey data
#' @param control.data data simulated from control distribution
#' @param case.data data simulated from case distribution
#' @param type method; likelihood, youden or specificity
#' @param nboot number of replicates in bootstrap
#'
#' @return
#' @export
#'
#' @examples
bootstrap.analyze <- function( data.serosurvey, control.data, case.data , type = 'likelihood', nboot = 1000){

  result.boot <- rep(0, nboot)

  for (i in seq(1, nboot, 1)){
    data.serosurvey.resample <- sample(x = data.serosurvey, size = length(data.serosurvey), replace = TRUE)
    case.data.resample <- sample(x = case.data, size = length(case.data), replace = TRUE)
    control.data.resample <- sample(x = control.data, size = length(control.data), replace = TRUE)

    if (type == 'likelihood'){
      pdfs <- sim.validation.data(control.data = control.data.resample, case.data = case.data.resample)
      result.boot[i] <- analyze.serosurvey.emperical(data.serosurvey = data.serosurvey.resample, pdf.case.smooth = pdfs[[1]],
                                                     pdf.control.smooth = pdfs[[2]], type = 'likelihood')[[1]]
    } else if (type == 'MY'){
      max.youden.emp <- observed.threshold(control.data.resample, case.data.resample, method = 'youden')
      output.clas <- analyze.serosurvey.emperical(data.serosurvey = data.serosurvey.resample, type = 'youden', threshold = max.youden.emp[1],
                                            TNR = max.youden.emp[3], TPR = max.youden.emp[2], correct = T)[1]
      result.boot[i] <-  output.clas[[1]]
    } else if (type == 'HS'){
      high.sens.emp <- observed.threshold(control.data.resample, case.data.resample, method = 'specificity')
      output.spec <- analyze.serosurvey.emperical(data.serosurvey = data.serosurvey.resample, type = 'youden', threshold = high.sens.emp[1],
                                            TNR = high.sens.emp[3], TPR = high.sens.emp[2], correct = T)[1]
      result.boot[i] <- output.spec[[1]]
    }

  }

  result.boot <- sort(result.boot)

  mean <- mean(result.boot)
  lower.CI <- result.boot[floor(nboot*2.5/100)]
  upper.CI <- result.boot[ceiling(nboot*97.5/100)]

  return(list(mean, lower.CI, upper.CI))

}

#' simulate.serosurvey.MS
#'
#' @param N_tests number of tests in the serosurvey
#' @param mean.S shape of severe case distribution
#' @param mean.M shape of mild case distribution
#' @param prev true cumulative incidence
#' @param p.M proportion of cases that is mild
#'
#' @return
#' @export
#'
#' @examples
simulate.serosurvey.MS <- function( N_tests = 1e4, mean.S = 6, mean.M = 4, prev = 0.08, p.M = 0.86 ){

  # determine how many of the tests are positive
  N_pos <- rbinom(N_tests,1,prev)

  titer.data <- rep(0, N_tests)

  # create matrix with titer values of all positive and negative values
  for (j in seq(1,N_tests,1)){
    if (N_pos[j] == 1){
      if (rbinom(1,1,p.M)){
        titer.data[j] <- rgamma(1, shape = mean.M, scale = 1)
      } else {
        titer.data[j] <- rgamma(1, shape = mean.S, scale = 1)
      }
    } else if (N_pos[j] == 0){
      titer.data[j] <- rgamma(1, shape = 1, scale = 1)
    }
  }

  return(titer.data)
}

#' analyze.serosurvey.AS.emp
#'
#' @param data.serosurvey simulated serosurvey data
#' @param prev.mild proportion of cases that is mild
#' @param pdf.control.smooth estimated distribution of control sera
#' @param pdf.case.smooth estimated distribution of case sera
#' @param type method; likelihood, youden or specificity
#' @param mean.M shape of mild case distribution
#' @param threshold cutoff-value for cutoff-based measures
#' @param TNR true negative rate
#' @param TPR true positive rate
#' @param correct whether the ad-hoc Rogan-Gladen should be applied
#' @param complex  whether the ad-hoc Rogan-Gladen should be applied
#' @param case.data data of case sera
#'
#' @return
#' @export
#'
#' @examples
analyze.serosurvey.AS.emp <- function(data.serosurvey = data.serosurvey, prev.mild = 0.86*0.1, pdf.control.smooth, pdf.case.smooth,
         type = 'likelihood', mean.M = NULL, threshold = NULL, TNR = NULL, TPR = NULL, correct = FALSE, complex = TRUE, case.data = case.data){

  if (type == 'likelihood'){
    # initial guess of prevalence
    prev_init <- 0.1 #sum(data.serosurvey>(max(data.serosurvey)/3))/length(data.serosurvey)

    if (complex){
      fit <- try(optim(par = c(prev.tot = prev_init, prev.mild = 0.8*prev_init, mean.M =8), fn = function(x){ ll.titer.MS.emp(prev.tot=x[1], prev.mild=x[2], mean.M=x[3],
                                                                                                                              pdf.case.smooth = pdf.case.smooth,
                                                                                          pdf.control.smooth = pdf.control.smooth, data.serosurvey = data.serosurvey, both=TRUE, case.data=case.data) } ,
                       control = list(fnscale = -1), hessian = TRUE, lower = c(0.005,0.005,4), upper = c(1,1,12), method = "L-BFGS-B" )) # try MLE2 instead of optim
    } else {
      fit <- try(optim(par = c(prev.tot = prev_init), fn = function(x){ ll.titer.MS.emp(prev.tot=x[1], prev.mild=NULL, mean.M=NULL, pdf.case.smooth = pdf.case.smooth,
                                                                                        pdf.control.smooth = pdf.control.smooth,  data.serosurvey = data.serosurvey, both=FALSE) } ,
                       control = list(fnscale = -1), method = 'Brent', lower = 0, upper = 1, hessian = TRUE)) # try MLE2 instead of optim
    }

    prev_est <- fit$par
    return(c(prev_est, ll=fit$value))
  } else if ((type == 'naive') & correct){
    prev_est <- (sum(data.serosurvey>threshold)/length(data.serosurvey) + TNR - 1  ) / ( TPR + TNR - 1 )
    return(c(prev_est))
  } else if ((type == 'naive') & !correct){
    prev_est <- (sum(data.serosurvey>threshold)/length(data.serosurvey) )
    return(c(prev_est))
  }

}

#' ll.titer.MS.emp
#'
#' @param prev.tot total cumulative incidence
#' @param prev.mild cumulative incidence
#' @param mean.M shape of mild case distribution
#' @param pdf.control.smooth estimated
#' @param pdf.case.smooth estimated distribution of case sera
#' @param data.serosurvey simulated serosurvey data
#' @param both whether both distributions should be included
#' @param case.data data simulated from case distribution
#'
#' @return
#' @export
#'
#' @examples
ll.titer.MS.emp <- function(prev.tot, prev.mild, mean.M, pdf.control.smooth, pdf.case.smooth, data.serosurvey, both, case.data){

  # function calculating the likelihood of observing this data given the prevalence
  if (both){
    if (prev.mild < prev.tot){
      l <- prev.mild * dgamma(data.serosurvey, shape = mean.M, scale = 1) + #pdf.asym.smooth( data.serosurvey) +
        (prev.tot-prev.mild) *pdf.case.smooth( data.serosurvey )+  #pdf.case.smooth( data.serosurvey )
        (1-prev.tot) *pdf.control.smooth( data.serosurvey )# dgamma(data.serosurvey, shape = 1, scale = 1) ## #

      ll <- sum(log(l)) # check log1p function
    } else {
      ll <- -1e15
    }
  } else {
    l <-   prev.tot * pdf.case.smooth(data.serosurvey) +
      (1-prev.tot) * pdf.control.smooth( data.serosurvey )
    ll <- sum(log(l)) # check log1p function

  }

  return(ll)
}
