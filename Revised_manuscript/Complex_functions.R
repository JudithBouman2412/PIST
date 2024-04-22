#' RunBayesian_models
#'
#' @param N number of replicates
#' @param true_prevs true cumulative incidences modeled
#' @param means means of the case distributions
#' @param AUCs AUC-ROC values of the case and control distribution
#' @param file.save file where to save the data
#' @param rstan.model the rstan model to run
#'
#' @return
#' @export
#'
#' @examples
RunBayesian_models <- function(N=50, true_prevs=c(0.01,0.04,0.08), means = c(2.322645, 2.719439, 3.314629, 4.306613, 5.320641, 6.643287),
                               AUCs = c(0.801,0.849,0.9,0.95,0.975,0.990),
                               file.save = 'data.bayesian.Fig1.txt', rstan.model = "/Users/judith/testing/Revised_manuscript/cutoff_based.stan"){
  #create matrices for results
  final <- matrix(NA, nrow = 2, ncol = 7)
  final[,4] <- c("Max Youden - Bayesian", "High specificity - Bayesian")

  M_cutoff_based = stan_model(file=rstan.model)

  for (p in seq(1, length(true_prevs), 1)){
    final[,3] <- rep(true_prevs[p], 2)

    #simulate results for likelihood method
    for (i in seq(1, length(means), 1)){
      final[,2] <- rep( round(AUCs[i], 3), 2)

      for (j in seq(1,N,1)){
        final[,3] <- rep(true_prevs[p], 2)
        # Simulate validation data
        control.data <- rgamma(5000, shape = 1, scale = 1)
        case.data <- rgamma(5000, shape = means[i], scale = 1)

        max.youden.emp <- observed.threshold(control.data, case.data, method = 'youden')
        high.sens.emp <- observed.threshold(control.data, case.data, method = 'specificity')

        # simulate data given those distributions
        data.serosurvey <- sim.serosurvey( N_tests = N_tests,  mean = means[i], prev=true_prevs[p])

        # Create dichotomized data
        case.data.dich.clas <- case.data
        case.data.dich.clas[case.data<max.youden.emp[1]] <- 0
        case.data.dich.clas[case.data>=max.youden.emp[1]] <- 1

        control.data.dich.clas <- control.data
        control.data.dich.clas[control.data<max.youden.emp[1]] <- 0
        control.data.dich.clas[control.data>=max.youden.emp[1]] <- 1

        data.serosurvey.dich.clas <- data.serosurvey
        data.serosurvey.dich.clas[data.serosurvey.dich.clas<max.youden.emp[1]] <- 0
        data.serosurvey.dich.clas[data.serosurvey.dich.clas>=max.youden.emp[1]] <- 1

        # create dataset for Bayesian model
        sim_list_clas = list(N1=length(case.data.dich.clas),
                             N2=length(control.data.dich.clas),
                             N3=length(data.serosurvey.dich.clas),
                             cutoff_case=case.data.dich.clas,
                             cutoff_control=control.data.dich.clas,
                             cutoff_serosurvey=data.serosurvey.dich.clas,
                             prior_sensitivity=c(1,1), # gamma prior for shape parameters
                             prior_specificity=c(1,1)) # gamma prior for rate parameters

        # with bayesian correction from Larremore/Carpenter
        S_cutoff_based_clas = sampling(M_cutoff_based,data=sim_list_clas,chains=4,iter=2000)
        bh_summary_clas <- summary(S_cutoff_based_clas)$summary %>%
          as_tibble() %>%
          mutate(variable = rownames(.)) %>%
          select(variable, everything()) %>%
          as_data_frame()

        final[1,1] <- bh_summary_clas$mean[3]
        final[1,5] <- bh_summary_clas$mean[3]
        final[1,6] <- bh_summary_clas$`2.5%`[3]
        final[1,7] <- bh_summary_clas$`97.5%`[3]

        # Create dichotomized data
        case.data.dich.spec <- case.data
        case.data.dich.spec[case.data<high.sens.emp[1]] <- 0
        case.data.dich.spec[case.data>=high.sens.emp[1]] <- 1

        control.data.dich.spec <- control.data
        control.data.dich.spec[control.data<high.sens.emp[1]] <- 0
        control.data.dich.spec[control.data>=high.sens.emp[1]] <- 1

        data.serosurvey.dich.spec <- data.serosurvey
        data.serosurvey.dich.spec[data.serosurvey.dich.spec<high.sens.emp[1]] <- 0
        data.serosurvey.dich.spec[data.serosurvey.dich.spec>=high.sens.emp[1]] <- 1

        # create dataset for Bayesian model
        sim_list_spec = list(N1=length(case.data.dich.spec),
                             N2=length(control.data.dich.spec),
                             N3=length(data.serosurvey.dich.spec),
                             cutoff_case=case.data.dich.spec,
                             cutoff_control=control.data.dich.spec,
                             cutoff_serosurvey=data.serosurvey.dich.spec,
                             prior_sensitivity=c(1,1), # gamma prior for shape parameters
                             prior_specificity=c(1,1)) # gamma prior for rate parameters

        # with bayesian correction from Larremore/Carpenter
        S_cutoff_based_spec = sampling(M_cutoff_based, data=sim_list_spec,chains=4,iter=2000)

        bh_summary_spec <- summary(S_cutoff_based_spec)$summary %>%
          as_tibble() %>%
          mutate(variable = rownames(.)) %>%
          select(variable, everything()) %>%
          as_data_frame()

        final[2,1] <- bh_summary_spec$mean[3]
        final[2,5] <- bh_summary_spec$mean[3]
        final[2,6] <- bh_summary_spec$`2.5%`[3]
        final[2,7] <- bh_summary_spec$`97.5%`[3]

        write.table(final, file = file.save, append = TRUE, col.names = FALSE )
      }

    }

  }
}


#' RunNonBayesian_models
#'
#' @param N number of replicates
#' @param true_prevs true cumulative incidences modeled
#' @param means means of the case distributions
#' @param AUCs AUC-ROC values of the case and control distribution
#' @param file.save file where to save the data
#'
#' @return
#' @export
#'
#' @examples
RunNonBayesian_models <- function(N=50, true_prevs=c(0.01,0.04,0.08), means = c(2.322645, 2.719439, 3.314629, 4.306613, 5.320641, 6.643287),
                                  AUCs = c(0.801,0.849,0.9,0.95,0.975,0.990), file.save = 'data.nonbayesian.Fig1.txt'){

  final <- matrix(NA, nrow = 3, ncol = 7)
  N_tests = 1e4

  for (p in seq(1, length(true_prevs), 1)){
    final[,3] <- rep(true_prevs[p], 3)

    #simulate results for likelihood method
    for (i in seq(1, length(means), 1)){
      final[,2] <- rep( round(AUCs[i], 3), 3)

      for (j in seq(1,N,1)){
          final[,4] <- c("Mixture model" , "Max Youden - RG", "High specificity - RG")

          # Simulate validation data
          control.data <- rgamma(5000, shape = 1, scale = 1)
          case.data <- rgamma(5000, shape = means[i], scale = 1)

          # simulate data given those distributions
          data.serosurvey <- sim.serosurvey( N_tests = N_tests,  mean =  means[i], prev=true_prevs[p])

          ############### likelihood-based method #######################
          # simulate validation data and create pdf for both distributions from there
          pdfs <- sim.validation.data(control.data = control.data, case.data = case.data)

          # analyze data and estimate prev
          output <- analyze.serosurvey.emperical(data.serosurvey = data.serosurvey, pdf.case.smooth = pdfs[[1]],
                                                 pdf.control.smooth = pdfs[[2]], type = 'likelihood')

          # what is the probability that our estimated prevalence is equal to the true prevalence?
          final[1,1] <- output[[1]]

          boot.output <- bootstrap.analyze( data.serosurvey=data.serosurvey, control.data=control.data,
                                            case.data=case.data , type = 'likelihood', nboot = 1000)
          final[1,5] <- boot.output[[1]]
          final[1,6] <- boot.output[[2]]
          final[1,7] <- boot.output[[3]]

          ############### cutoff-based methods #######################

          # calculate thresholds based on observed validation data
          max.youden.emp <- observed.threshold(control.data, case.data, method = 'youden')
          high.sens.emp <- observed.threshold(control.data, case.data, method = 'specificity')

          # analyze data and estimate prev, with classical correction
          output.clas <- analyze.serosurvey.emperical(data.serosurvey = data.serosurvey, type = 'youden', threshold = max.youden.emp[1],
                                                TNR = max.youden.emp[3], TPR = max.youden.emp[2], correct = T)[1]

          output.spec <- analyze.serosurvey.emperical(data.serosurvey = data.serosurvey, type = 'youden', threshold = high.sens.emp[1],
                                                TNR = high.sens.emp[3], TPR = high.sens.emp[2], correct = T)[1]

          # what is the probability that our estimated prevalence is equal to the true prevalence?
          final[2,1] <- output.clas[[1]]

          boot.clas <- bootstrap.analyze( data.serosurvey=data.serosurvey, control.data=control.data,
                                          case.data=case.data , type = 'MY', nboot = 1000)
          final[2,5] <- boot.clas[[1]]
          final[2,6] <- boot.clas[[2]]
          final[2,7] <- boot.clas[[3]]

          boot.spec <- bootstrap.analyze( data.serosurvey=data.serosurvey, control.data=control.data,
                                          case.data=case.data , type = 'HS', nboot = 1000)
          final[3,5] <- boot.spec[[1]]
          final[3,6] <- boot.spec[[2]]
          final[3,7] <- boot.spec[[3]]

          final[3,1] <- output.spec[[1]]

          #save final to result file
          write.table(final, file = file.save, append = TRUE, col.names = FALSE )
      }
    }
  }

}

#' power
#'
#' @param power.data
#' @param sample.sizes sample size of the cases and controls in the validation data
#' @param AUC tested AUC values
#' @param serosurvey.size size of the serosurvey
#'
#' @return
#' @export
#'
#' @examples
power <- function(power.data, sample.sizes, AUC, serosurvey.size){
  # calculates the statistical power per data unity.

  # first of 10000 in the serosurvey
  power <- matrix(0, nrow = length(sample.sizes), ncol = length(AUC) )

  for (i in seq(1, length(sample.sizes),1)){
    for (j in seq(1, length(AUC),1)){
      data.now <- power.data[ ( power.data[,1]==AUC[j] & power.data[,4]==sample.sizes[i] & power.data[,2]==serosurvey.size),  ]
      #print(data.now)
      lower <-  data.now[1,3] - data.now[1,3]*0.25
      upper <- data.now[1,3] + data.now[1,3]*0.25
      tot <- dim(data.now)[1]
      print(tot)
      power[i,j] <- sum( (data.now[,5]>lower & data.now[,5]<upper) & ( data.now[,3]>(data.now[,6]-2*data.now[,6]) & data.now[,3]<(data.now[,5]+2*data.now[,6]))  )/tot
    }
  }
  return(power)
}
