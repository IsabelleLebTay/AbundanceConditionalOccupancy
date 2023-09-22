
data {
    int<lower=0> I; 
  int<lower=0, upper=1> y[I];  // Collapsed observed occurrence for each site
  
  // site-level covariates for occupancy
  vector[I] latitude;
  vector[I] longitude;
  vector[I] percent_conifer;
  vector[I] percent_deciduous;
  vector[I] percent_pine;
}

parameters {
  real beta_0;  // Intercept for occupancy probability
  real beta_latitude;
  real beta_longitude;
  real beta_percent_conifer;
  real beta_percent_deciduous;
  real beta_percent_pine;
}

model {
  // Priors for occupancy coefficients
  beta_0 ~ normal(0, 5);
  beta_latitude ~ normal(0, 5);
  beta_longitude ~ normal(0, 5);
  beta_percent_conifer ~ normal(0, 5);
  beta_percent_deciduous ~ normal(0, 5);
  beta_percent_pine ~ normal(0, 5);

  for (i in 1:I) {
    // Linear predictor for occupancy probability for site i
    real lin_pred_occupancy = beta_0 + 
                              beta_latitude * latitude[i] + 
                              beta_longitude * longitude[i];

    
    // Compute the occupancy probability for site i
    real psi_i = inv_logit(lin_pred_occupancy);
    
    // Bernoulli likelihood for site i
    y[i] ~ bernoulli(psi_i);
  }
}
