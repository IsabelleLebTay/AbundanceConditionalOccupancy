data{
    int<lower=0> I; // number of sites
    array[I] int<lower=0> M;
    vector[I] age;
    vector[I] size;
    vector[I] percent_conifer;
    vector[I] percent_pine;
    vector[I] percent_deciduous;

}

parameters{
    real alpha;
    real beta_age;
    real beta_size;
    real beta_conifer;
    real beta_pine;
    real beta_deciduous;
   
}


model{

    // likelihood
    M ~ poisson_log(alpha + beta_age*age 
                      + beta_size*size 
                      + beta_deciduous * percent_deciduous 
                      + beta_pine * percent_pine
                      + beta_conifer * percent_conifer);


    //priors
    alpha ~ normal(0, 0.5);
    beta_age ~ normal(0,1);
    beta_size ~ normal(0,1);
    beta_conifer ~ normal(0,1);
    beta_pine ~ normal(0,1);
    beta_deciduous ~ normal(0,1);  
}


generated quantities {
  int M_pred[I];  // predicted counts
  real log_lik[I]; // log-likelihood for each predicted count
  for (i in 1:I) {
    // Generate predicted count for each observation
    M_pred[i] = poisson_log_rng(alpha 
                                  + beta_age * age[i]
                                  + beta_deciduous * percent_deciduous[i]
                                  + beta_pine * percent_pine[i]
                                  + beta_deciduous * percent_deciduous[i]);

    log_lik[i] = poisson_log_lpmf(M[i] | alpha 
                                          + beta_age * age[i] 
                                          + beta_size * size[i] 
                                          + beta_deciduous * percent_deciduous[i]
                                          + beta_pine * percent_pine[i]
                                          + beta_deciduous * percent_deciduous[i]);
  }
}