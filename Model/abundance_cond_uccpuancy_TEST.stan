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
  vector[I] size;
  vector[I] age;

  // visit-level covariates
    matrix[I, J] time_of_day;  // Adjusted to be a matrix
    matrix[I, J] Julian_date;  // Adjusted to be a matrix
}

transformed data {
  int y[I, J];  // observed presence/absence, derived from the count N
  int M[I];       // mean bird count for each site

  for (i in 1:I) {

    M[i] = sum(N[i,]) / J;  // compute the mean count for site i

    for (j in 1:J) {
      y[i, j] = (N[i, j] > 0);
    }
  }
}

parameters {
  // Occupancy
  real beta0_psi;  // intercept for occupancy
  real beta1_psi;  // coefficient for latitude
  real beta2_psi;  // coefficient for longitude
  real beta_conifer_psi;
  real beta_deciduous_psi;
  real beta_pine_psi;
  real logit_psi[I];
  real<lower=0> sigma_psi;
  real<lower=0> sigma_p;
  real<lower=0> sigma_lambda;
  matrix[I, J] logit_p;  // log odds for each site i and visit j
  
  // Detection
  real beta0_p;    // intercept for detection probability
  real beta1_p;    // coefficient for time of day
  real beta2_p;    // coefficient for Julian date
  
  // Abundance
  real beta0_lambda;     // intercept for abundance
  real beta_size;  // coefficient for size
  real beta_age;   // coefficient for age
  real beta_conifer_lambda;
  real beta_deciduous_lambda;
  real beta_pine_lambda;
  real log_lambda[I];
}

transformed parameters {
    real psi[I];
    for (i in 1:I) {
        psi[i] = 1 / (1 + exp(-logit_psi[i]));

    }
}

model {
  // logit_psi: estimated parameter of the latent unobserved ocucpancy
  for (i in 1:I) {
    logit_psi[i] ~ beta(beta0_psi + beta1_psi * latitude[i] + beta2_psi * longitude[i] + 
														beta_conifer_psi * percent_conifer[i] + beta_deciduous_psi * percent_deciduous[i] + 
														beta_pine_psi * percent_pine[i], sigma_psi);
		for (j in 1:J) {
															
            logit_p[i, j] ~ normal(beta0_p + beta1_p * time_of_day[i, j] + beta2_p * Julian_date[i, j], sigma_p);
        }
			
      // this should be a deterministic equation
      log_lambda[i] ~ normal(beta0_lambda + beta_size * size[i] + beta_age * age[i] + 
														 beta_conifer_lambda * percent_conifer[i] + beta_deciduous_lambda * percent_deciduous[i] + 
														 beta_pine_lambda * percent_pine[i], sigma_lambda);
			
		
  }

  for (i in 1:I) {
    //z[i] ~ bernoulli_logit(logit_psi[i]);
    target += log_mix(1 / (1 + exp(-logit_psi[i])), 0, 1);  // This is a trick to marginalize over the binary z

    for (j in 1:J) {
      y[i, j] ~ bernoulli_logit( psi[i] * logit_p[i, j]);
   
    }
    M[i] ~ poisson_log(log_lambda[i] + log(psi[i]));
  }
  
  // Priors

  // Detection (p) 
  beta0_p ~ normal(0, 10);
  beta1_p ~ normal(0, 10);
  beta2_p ~ normal(0, 10);

  // Occupancy (psi)
  beta0_psi ~ normal(0, 10);
  beta1_psi ~ normal(0, 10);
  beta2_psi ~ normal(0, 10);
  beta_conifer_psi ~ normal(0, 10);
  beta_deciduous_psi ~ normal(0, 10);
  beta_pine_psi ~ normal(0, 10);

  // Abundance (lambda)
  beta0_lambda  ~ normal(0, 10);
  beta_size ~ normal(0, 10);
  beta_age ~ normal(0, 10);

  beta_conifer_lambda ~ normal(0, 10);
  beta_deciduous_lambda ~ normal(0, 10);
  beta_pine_lambda ~ normal(0, 10);

  sigma_psi ~ normal(0, 2.5);
  sigma_p ~ normal(0, 2.5);
  sigma_lambda ~ normal(0, 2.5);
}


generated quantities {
  // declare parameters
  real logit_psi_pred[I];
  real logit_p_pred[I, J]; 
  real log_lambda_pred[I];
  int z_pred[I];
  int M_pred[I];
  int y_rep[I, J]; // posterior predictive of observed presence/absence


	for (i in 1:I) {
		for (j in 1:J) {
			logit_psi_pred[i] = beta0_psi + beta1_psi * latitude[i] + beta2_psi * longitude[i] + 
														beta_conifer_psi * percent_conifer[i] + beta_deciduous_psi * percent_deciduous[i] + 
														beta_pine_psi * percent_pine[i];
														
      logit_p_pred[i, j] = beta0_p + beta1_p * time_of_day[i, j] + beta2_p * Julian_date[i, j];
			
		}

        log_lambda_pred[i] = beta0_lambda + beta_size * size[i] + beta_age * age[i] + 
														 beta_conifer_lambda * percent_conifer[i] + beta_deciduous_lambda * percent_deciduous[i] + 
														 beta_pine_lambda * percent_pine[i];

    
	}
   for (i in 1:I) {
    z_pred[i] = bernoulli_logit_rng(logit_psi_pred[i]);

    for (j in 1:J) {
    // Generate a new presence/absence data point
  y_rep[i, j] = bernoulli_logit_rng( z_pred[i] * logit_p_pred[i, j]);

  M_pred[i] = poisson_log_rng(log_lambda_pred[i] * z_pred[i]);

    }
    }

}