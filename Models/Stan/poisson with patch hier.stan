data {
  int<lower=0> I; // number of sites
  int<lower=1, upper=2> patch[I]; // Patch indicator (1 or 2 instead of 0 or 1)
  vector[I] age;           
  array[I] int<lower=0> M;
  vector[I] size;
  vector[I] percent_conifer;
  vector[I] percent_pine;
  vector[I] percent_deciduous;
}

parameters {
  vector[2] alpha;           // Intercepts for each patch
  real<lower=0> sigma_alpha; // Standard deviation of intercepts
  real beta_age;           
  real beta_size;
  real beta_conifer;
  real beta_pine;
  real beta_deciduous;
  real beta_age_size;
}

model {
  // Priors
  alpha ~ normal(0, 1);
  sigma_alpha ~ exponential(1);
  beta_age ~ normal(0, 1);
  beta_size ~ normal(0,1);
  beta_age_size ~ normal(0,1);
  beta_conifer ~ normal(0,1);
  beta_pine ~ normal(0,1);
  beta_deciduous ~ normal(0,1);  

  // Likelihood

  M ~ poisson_log(alpha[patch] + beta_age * age
                      + beta_age_size .* age // .* indicates element-wise multiplication
                      + beta_size * size 
                      + beta_deciduous * percent_deciduous 
                      + beta_pine * percent_pine
                      + beta_conifer * percent_conifer);

//  for (i in 1:I) {
//    M[i] ~ poisson_log(alpha[patch[i]] + beta_age * age[i]);
//  }
}



/* data{
    int<lower=0> I; // number of sites
    array[I] int<lower=0> M;
    vector[I] age;
    vector[I] size;
    vector[I] percent_conifer;
    vector[I] percent_pine;
    vector[I] percent_deciduous;
    array[I] int<lower=0, upper = 1> patch;
}

parameters{
    // real alpha;
    real beta_age;
    real beta_size;
    real beta_conifer;
    real beta_pine;
    real beta_deciduous;
    real beta_age_size;
    real beta_patch;
    // vector[2] intercept;
    real<lower=0> sigma;
    real intercept[2;]
}


model{
    // likelihood
  
    M ~ poisson_log(intercept[patch] + beta_age * age 
                      + beta_age_size .* age // .* indicates element-wise multiplication
                      + beta_size * size 
                      + beta_deciduous * percent_deciduous 
                      + beta_pine * percent_pine
                      + beta_conifer * percent_conifer);
  
    
    
    
      // Priors
  alpha ~ normal(0, 10);
  sigma_alpha ~ cauchy(0, 2.5);
  beta_age ~ normal(0, 10);

  // Likelihood
  for (i in 1:N) {
    M[i] ~ poisson_log(alpha[patch[i]] + beta_age * age[i]);
    
    
    // intercept ~ normal(mu[patch], sigma);

    //priors
    // mu ~ uniform(1, 1);

    //sigma ~ exponential(1);
    
    beta_age ~ normal(0,1);
    beta_size ~ normal(0,1);
    beta_age_size ~ normal(0,1);
    beta_conifer ~ normal(0,1);
    beta_pine ~ normal(0,1);
    beta_deciduous ~ normal(0,1);  
    
}


 generated quantities {
  array[I] int<lower=0> M_pred;
  real log_lik[I]; // log-likelihood for each predicted count
  for (i in 1:I) {
    // Generate predicted count for each observation
    M_pred[i] = poisson_log_rng(intercept[patch[i]]
                                  + beta_age * age[i]
                                  + beta_deciduous * percent_deciduous[i]
                                  + beta_pine * percent_pine[i]
                                  + beta_deciduous * percent_deciduous[i]);

    log_lik[i] = poisson_log_lpmf(M[i] | intercept[patch[i]] 
                                          + beta_age * age[i] 
                                          + beta_size * size[i] 
                                          + beta_deciduous * percent_deciduous[i]
                                          + beta_pine * percent_pine[i]
                                          + beta_deciduous * percent_deciduous[i]);
  }
} */