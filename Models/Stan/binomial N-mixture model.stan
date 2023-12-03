// Binomial mixture model with covariates
data {
    int<lower=0> I;       // Number of sites
    int<lower=0> J;       // Number of temporal replications
    int<lower=0> counts[I, J]; // Counts, for each site and visit (observed)
    int<lower=0> K;       // Upper bound of population size, total study

    matrix[I, J] tod;  // Time of day for each site and visit
    matrix[I, J] toy;  // Julian date for each site and visit
    vector[I] age;
    vector[I] size;
    vector[I] percent_conifer;
    vector[I] percent_pine;
    vector[I] percent_deciduous;
}

transformed data {
  int<lower=0> max_count[I]; // per site

  for (i in 1:I)
    max_count[i] = max(counts[i]);
}

parameters {
    real alpha_0;  // Intercept for detection probability
    real alpha_1;  // Coefficient for time_of_day
    real alpha_2;  // Coefficient for Julian_date

    real beta0;
    real beta_age;
    real beta_size;
    real beta_conifer;
    real beta_pine;
    real beta_deciduous;
    real beta_age_size;
}

transformed parameters {
  vector[I] log_lambda; // Log population size
  // matrix[I, J] logit_p; // Logit detection probability

  log_lambda = beta0 + beta_age * age 
                      + beta_age_size .* age // .* indicates element-wise multiplication
                      + beta_size * size 
                      + beta_deciduous * percent_deciduous 
                      + beta_pine * percent_pine
                      + beta_conifer * percent_conifer;

  // logit_p = rep_matrix(alpha_0 + alpha_1 * tod + alpha_2 + toy, J);
}

model {
      // Priors for detection probability coefficients
    alpha_0 ~ normal(0, 1);
    alpha_1 ~ normal(0, 1);
    alpha_2 ~ normal(0, 1);

    //priors, abundance
    beta0 ~ normal(0, 0.5);
    beta_age ~ normal(0,1);
    beta_size ~ normal(0,1);
    beta_age_size ~ normal(0,1);
    beta_conifer ~ normal(0,1);
    beta_pine ~ normal(0,1);
    beta_deciduous ~ normal(0,1);  

  // Likelihood
  for (i in 1:I) {
    vector[K - max_count[i] + 1] log_prob_f; // lp is the log probability function. The vector of size Total_max - local_max +1

    for (n in 1:(K - max_count[i] + 1)) {
      log_prob_f[n] = poisson_log_lpmf(max_count[i] + n - 1 | log_lambda[i])
             + binomial_logit_lpmf(counts[i] | max_count[i] + n - 1, logit_p[i])}; // logit_p is the detection probability at that site
    target += log_sum_exp(log_prob_f);
  }
}

generated quantities {
    // I want to track lambda!
  int N[I];
  int totalN;

  for (i in 1:R)
    N[i] = poisson_log_rng(log_lambda[i]); 
    log_lik[i] = poisson_log_lpmf(N[i] | log_lambda[i]);
  totalN = sum(N);
}