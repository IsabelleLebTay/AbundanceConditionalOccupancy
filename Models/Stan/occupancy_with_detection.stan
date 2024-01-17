

data {
  int<lower=0> I;  // number of sites
  int<lower=0> J;  // number of visits per site
  matrix[I, J] N;     // N is the bird count for each site and visit (observed)
  matrix[I, J] time_of_day;  // Time of day for each site and visit
  matrix[I, J] Julian_date;  // Julian date for each site and visit
  
  // site-level covariates for occupancy
  vector[I] latitude;
  vector[I] longitude;
/*  vector[I] percent_conifer;
  vector[I] percent_deciduous;
  vector[I] percent_pine; */
  vector[I] age;
}

transformed data {
  int<lower=0, upper=1> y[I, J];  // matrix of size I, J of the observed presence/absence, derived from the count N
  for (i in 1:I) {
    for (j in 1:J) {
      y[i, j] = (N[i, j] > 0);
    }
  }
}

parameters {
  real beta_0;  // Intercept for occupancy probability
  real beta_latitude;
  real beta_longitude;
  
 /* real beta_percent_conifer;
  real beta_percent_deciduous;
  real beta_percent_pine; */
  
  real beta_age; // Coefficient for the age

  real alpha_0;  // Intercept for detection probability
  real alpha_1;  // Coefficient for time_of_day
  real alpha_2;  // Coefficient for Julian_date
}

model {
  // Priors for occupancy coefficients
  beta_0 ~ normal(0, 1);
  beta_latitude ~ normal(0, 1);
  beta_longitude ~ normal(0, 1);
  beta_age ~ normal(0,1);

 /* beta_percent_conifer ~ normal(0, 1);
  beta_percent_deciduous ~ normal(0, 1);
  beta_percent_pine ~ normal(0, 1); */
  

  // Priors for detection probability coefficients
  alpha_0 ~ normal(0, 1);
  alpha_1 ~ normal(0, 1);
  alpha_2 ~ normal(0, 1);

  for (i in 1:I) {
    // Linear predictor for occupancy probability for site i

    real lin_pred_occupancy = beta_0; // intercept only occupancy model

    
    // Compute the occupancy probability for site i
    real psi_i = inv_logit(lin_pred_occupancy);
    
    for (j in 1:J) {
      // Linear predictor for detection probability at visit j for site i
      real lin_pred_detection = alpha_0 
                                + alpha_1 * time_of_day[i,j] 
                                + alpha_2 * Julian_date[i,j]
                                + beta_latitude * latitude[i] 
                                + beta_longitude * longitude[i];
      
      // Compute the detection probability at visit j for site i
      real p_ij = inv_logit(lin_pred_detection);
      
      // Log-likelihood when z_i = 0 //(species not present at site i)
      real log_lik_0 = bernoulli_lpmf(y[i,j] | 0);
      
      // Log-likelihood when z_i = 1 (species present at site i)
      real log_lik_1 = bernoulli_lpmf(y[i,j] | p_ij);
      
      // Marginalize out z_i by adding the log-likelihoods for site i
      target += log_sum_exp(log_lik_0, log_lik_1 + log(psi_i));
    }
  }
}

