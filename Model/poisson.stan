data{
    int<lower=0> I; // number of sites
    array[I] int<lower=0> M;
    vector[I] age;
    vector[I] size;

}

parameters{
    real alpha;
    real beta_age;
    real beta_size;

   
}


model{

    // likelihood
    M ~ poisson_log(alpha + beta_age*age + beta_size*size);


    //priors
    alpha ~ normal(0, 0.5);
    beta_age ~ normal(0,1);
    beta_size ~ normal(0,1);

    
}


generated quantities {
  int M_pred[I];  // predicted counts
  real log_lik[I]; // log-likelihood for each predicted count
  for (i in 1:I) {
    // Generate predicted count for each observation
    M_pred[i] = poisson_log_rng(alpha + beta_age * age[i] + beta_size*size[i]);
    log_lik[i] = poisson_log_lpmf(M[i] | alpha + beta_age * age[i] + beta_size*size[i]);
  }
}