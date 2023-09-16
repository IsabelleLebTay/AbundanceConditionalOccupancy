data {
  int<lower=0> I;  // number of sites (414 total)
  int<lower=0> J;  // number of visits per site (15)
  int N[I, J];     // N is the bird count for each site and visit (observed)

  // covariates: observed continuous values
  // site-level covariates
  vector[I] latitude;
  vector[I] longitude;
  vector[I] percent_conifer;
  vector[I] percent_deciduous;
  vector[I] percent_pine;
  vector[I] percent_mixed;

  // visit-level covariates
  vector[J] time_of_day;
  vector[J] Julian_date;
  vector[J] temperature;
}

parameters {
  // Occupancy
  real beta0_psi;  // intercept for occupancy
  real beta1_psi;  // coefficient for latitude
  real beta2_psi;  // coefficient for longitude
  real beta_conifer_psi;
  real beta_deciduous_psi;
  real beta_pine_psi;
  real beta_mixed_psi;
  real z[I];
  
  // Detection
  real beta0_p;    // intercept for detection probability
  real beta1_p;    // coefficient for time of day
  real beta2_p;    // coefficient for Julian date
  real beta3_p;    // coefficeint for temperature
  
  // Abundance
  real beta0_lambda;     // intercept for abundance
  real beta_size;  // coefficient for size
  real beta_age;   // coefficient for age
  real beta_conifer_lambda;
  real beta_deciduous_lambda;
  real beta_pine_lambda;
  real beta_mixed_lambda;
  real M[I];
}

transformed data {
  int y[I, J];  // observed presence/absence, derived from the count N
  real M[I];       // mean bird count for each site

  for (i in 1:I) {

    M[i] = sum(N[i,]) / J;  // compute the mean count for site i
    for (j in 1:J) {
      y[i, j] = (N[i, j] > 0);
    }
  }
}

model {
  // logit_psi: estimated parameter of the latent unobserved ocucpancy
  vector[I] logit_psi = beta0_psi + beta1_psi * latitude + beta2_psi * longitude + 
                        beta_conifer_psi * percent_conifer + beta_deciduous_psi * percent_deciduous + 
                        beta_pine_psi * percent_pine + beta_mixed_psi * percent_mixed;
                        
  
  // logit_p: estimated detection probability for each site and visit
  matrix[I, J] logit_p = beta0_p + beta1_p * time_of_day + beta2_p * Julian_date + beta3_p * temperature;

  // log_lambda: estimated rate in the abundance function

  vector[I] log_lambda = beta0_lambda + beta_size * size + beta_age * age + 
                         beta_conifer_lambda * percent_conifer + beta_deciduous_lambda * percent_deciduous + 
                         beta_pine_lambda * percent_pine + beta_mixed_lambda * percent_mixed;

  for (i in 1:I) {
    z[i] ~ bernoulli_logit(logit_psi[i]);
    for (j in 1:J) {
      y[i, j] ~ bernoulli_logit( z[i] * logit_p[i, j]);
      
      M[i] ~ poisson_log(log_lambda[i] * z[i]);
      
    }

  }
  
  // Priors

  // Detection (p) 
  beta0_p ~ normal(0, 10);
  beta1_p ~ normal(0, 10);
  beta2_p ~ normal(0, 10);
  beta3_p ~ normal(0, 10);

  // Occupancy (psi)
  beta0_psi ~ normal(0, 10);
  beta1_psi ~ normal(0, 10);
  beta2_psi ~ normal(0, 10);
  beta_conifer_psi ~ normal(0, 10);
  beta_deciduous_psi ~ normal(0, 10);
  beta_pine_psi ~ normal(0, 10);
  beta_mixed_psi ~ normal(0, 10);

  // Abundance (lambda)
  beta0_lambda  ~ normal(0, 10);
  beta_size ~ normal(0, 10);
  beta_age ~ normal(0, 10);

  beta_conifer_lambda ~ normal(0, 10);
  beta_deciduous_lambda ~ normal(0, 10);
  beta_pine_lambda ~ normal(0, 10);
  beta_mixed_lambda ~ normal(0, 10);
}


generated quantities {
  // declare parameters
  real logit_psi_pred[I];
  real logit_p_pred[I, J]; 
  real log_lambda_pred[I];
  real z_pred[I];
  real M_pred[I];


	int y_rep[I, J]; // posterior predictive of observed presence/absence
	int N_rep[I, J]; // posterior predictive of observed counts

	for (i in 1:I) {
		for (j in 1:J) {
			logit_psi_pred[i] = beta0_psi + beta1_psi * latitude[i] + beta2_psi * longitude[i] + 
														beta_conifer_psi * percent_conifer[i] + beta_deciduous_psi * percent_deciduous[i] + 
														beta_pine_psi * percent_pine[i] + beta_mixed_psi * percent_mixed[i];
														
      logit_p_pred[i, j] = beta0_p + beta1_p * time_of_day[j] + beta2_p * Julian_date[j];
			
      log_lambda_pred[i] = beta0_lambda + beta_size * size[i] + beta_age * age[i] + 
														 beta_conifer_lambda * percent_conifer[i] + beta_deciduous_lambda * percent_deciduous[i] + 
														 beta_pine_lambda * percent_pine[i] + beta_mixed_lambda * percent_mixed[i];
			
		}
  
	}
   for (i in 1:I) {
    z_pred[i] ~ bernoulli_logit(logit_psi_pred[i]);
    for (j in 1:J) {
    // Generate a new presence/absence data point
  y_rep[i, j] = bernoulli_logit_rng( z_pred[i] * logit_p_pred[i, j]);

  M_pred[i] = poisson_log_rng(log_lambda_pred[i] * z_pred[i]);

    }
}

}