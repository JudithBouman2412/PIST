data {
  int N1; 
  int N2; 
  int N3; 
  int cutoff_case[N1];
  int cutoff_control[N2];
  int cutoff_serosurvey[N3];
  real prior_sensitivity[2];
  real prior_specificity[2];
}
parameters {
  real<lower=0,upper=1> sensitivity;
  real<lower=0,upper=1> specificity;
  simplex[2] prevalence;
}
model {  
  // prior distributions
  sensitivity ~ beta(prior_sensitivity[1],prior_sensitivity[2]);
  specificity ~ beta(prior_specificity[1],prior_specificity[2]);
  
  // likelihoods
  target += bernoulli_lpmf(cutoff_case|sensitivity);
  target += bernoulli_lpmf(cutoff_control|1-specificity);
  target += bernoulli_lpmf(cutoff_serosurvey|prevalence[1]*sensitivity+prevalence[2]*(1-specificity));
}
generated quantities {
  int pred_cutoff_case[N1];
  int pred_cutoff_control[N2];
  int pred_cutoff_serosurvey[N3];

  for(n in 1:N1) pred_cutoff_case[n] = bernoulli_rng(sensitivity); 
  for(n in 1:N2) pred_cutoff_control[n] = bernoulli_rng(1-specificity);
  for(n in 1:N3) pred_cutoff_serosurvey[n] = bernoulli_rng(prevalence[1]*sensitivity+prevalence[2]*(1-specificity));
}
