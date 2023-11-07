data {
  int<lower=0> I;  // number of sites (414 total)
  int<lower=0> J;  // number of visits per site (15)
  int N[I, J];     // N is the bird count for each site and visit (observed)
  int<lower=0> M[I]; // either rounded mean, or maximum, depending on what's being fed into the model.

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
  int<lower=0, upper=1> y[I, J];  // matrix of size I, J of the observed presence/absence, derived from the count N
  // int<lower=0> M[I];       // mean bird count for each site

  for (i in 1:I) {

    // M[i] = round(sum(N[i,]) / J);  // compute the mean count for site i, rounded so it works with the poisson_lpmf, Calcualting it in python first then feeding as data

    for (j in 1:J) {
      y[i, j] = (N[i, j] > 0);
    }
  }
}


parameters {
  // Occupancy
  real beta_0;  // Intercept for occupancy probability
  real beta_latitude;
  real beta_longitude;
  real beta_percent_conifer;
  real beta_percent_deciduous;
  real beta_percent_pine;

  // Detection
  real alpha_0;  // Intercept for detection probability
  real alpha_1;  // Coefficient for time_of_day
  real alpha_2;  // Coefficient for Julian_date

  // Abundance
  real gamma_0;  // Intercept for Poisson lambda
  real gamma_size;
  real gamma_age;
  real gamma_conifer_lambda;
  real gamma_deciduous_lambda;
  real gamma_pine_lambda;
}



model {
  for (i in 1:I) {
    // Linear predictor for occupancy probability for site i
    real lin_pred_occupancy = beta_0 + 
                              beta_latitude * latitude[i] + 
                              beta_longitude * longitude[i] + 
                              beta_percent_conifer * percent_conifer[i] + 
                              beta_percent_deciduous * percent_deciduous[i] +
                              beta_percent_pine * percent_pine[i];
    
    // Compute the occupancy probability for site i
    real psi_i = inv_logit(lin_pred_occupancy);
    
    // Linear predictor for Poisson lambda for site i
    real lin_pred_lambda = gamma_0 + gamma_size * size[i] + gamma_age * age[i] +
                                gamma_conifer_lambda * percent_conifer[i] + gamma_deciduous_lambda * percent_deciduous[i] + 
  														gamma_pine_lambda * percent_pine[i];
    
    // Compute Poisson lambda for site i
    real lambda_i = exp(lin_pred_lambda);
    
    // Zero-inflated Poisson likelihood for site i
    target += log_mix(psi_i, poisson_lpmf(M[i] | lambda_i), 0);
    
    for (j in 1:J) {
      // Linear predictor for detection probability at visit j for site i
      real lin_pred_detection = alpha_0 + alpha_1 * time_of_day[i,j] + alpha_2 * Julian_date[i,j];
      
      // Compute the detection probability at visit j for site i
      real p_ij = inv_logit(lin_pred_detection);
      
      // Log-likelihood when z_i = 0 (species not present at site i)
      real log_lik_0 = bernoulli_lpmf(y[i,j] | 0);
      
      // Log-likelihood when z_i = 1 (species present at site i)
      real log_lik_1 = bernoulli_lpmf(y[i,j] | p_ij);
      
      // Marginalize out z_i by adding the log-likelihoods for site i
      target += log_sum_exp(log_lik_0, log_lik_1 - log(psi_i));



  beta_0 ~ normal(0, 5);
  beta_latitude ~ normal(0, 5);
  beta_longitude ~ normal(0, 5);
  beta_percent_conifer ~ normal(0, 5);
  beta_percent_deciduous ~ normal(0, 5);
  beta_percent_pine ~ normal(0, 5);

  // Priors for Poisson lambda coefficients
  gamma_0 ~ normal(0, 5);
  gamma_size ~ normal(0, 5);
  gamma_age ~ normal(0, 5);

  // Priors for detection probability coefficients
  alpha_0 ~ normal(0, 5);
  alpha_1 ~ normal(0, 5);
  alpha_2 ~ normal(0, 5);
  gamma_pine_lambda ~ normal(0,5);
  gamma_conifer_lambda~ normal(0,5);
  gamma_deciduous_lambda~ normal(0,5);
    }
  }
}

/* model {
  //maths
  // z ~ bernoulli_logit(psi); //  estimated latent occupancy/ bernoulli_logit expects only one real argument
  // y ~ Bernoulli(z , p_detection); // vectorized, the likelihood of the observed data of presence/absence

  // not vectorized
  // for i in (1:I) {
  //  z[i] ~ bernoulli_logit(beta0_psi + beta1_psi * latitude[i] + beta2_psi * longitude[i] + 
	//													beta_conifer_psi * percent_conifer[i] + beta_deciduous_psi * percent_deciduous[i] + 
	//													beta_pine_psi * percent_pine[i], sigma_occupancy_process[i]);

  //  for j in (1:J) {
  //    p_detection[i, j] ~ bernoulli_logit(beta0_p + beta1_p * time_of_day[i, j] + beta2_p * Julian_date[i, j]); //
  //  }

  // vectorized likelihoods
  // occupancy
  // sampling of y, which is an array of ints of size J, vecotrized for each site (size I)
  y[J] ~ bernoulli_logit(z[I] * p_detection[J]);

  // to estimate z (the latent occupancy), we have to marginalise on the joint distribution of z so the HMC can sample from a gradient

  z ~ bernoulli_logit(beta0_psi + beta1_psi * latitude + beta2_psi * longitude + 
                        beta_conifer_psi * percent_conifer + beta_deciduous_psi * percent_deciduous + 
                        beta_pine_psi * percent_pine);
  for (j in 1:J) {
    p_detection[j] ~ beta(beta0_p + beta1_p * time_of_day[j] + beta2_p * Julian_date[j], sigma_ObservationProcess); //
  }



  // abundance, not vectorized
//    for i in (1:I) {
//      M[i] ~ poisson_log(z[i] * (beta0_lambda + beta_size * size[i] + beta_age * age[i] + 
//														 beta_conifer_lambda * percent_conifer[i] + beta_deciduous_lambda * percent_deciduous[i] + 
//														 beta_pine_lambda * percent_pine[i]));
//    }

  M ~ poisson_log(z * lambda);  // is this vectorized?
 
  lambda = poisson_log(beta0_lambda + beta_size * size + beta_age * age + 
                            beta_conifer_lambda * percent_conifer + beta_deciduous_lambda * percent_deciduous + 
                            beta_pine_lambda * percent_pine);


  // Priors
  // Detection (p) 
  beta0_p ~ normal(0, 1);
  beta1_p ~ normal(0, 10);
  beta2_p ~ normal(0, 10);
  sigma_ObservationProcess ~ cauchy(0,1);
  
  // Occupancy (psi)
  beta0_psi ~ normal(0, 1);
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
} */

/* generated quantities{
  // parameters
  int z_pred[I];
  vector[I] p_detection_pred[J]; 
  vector[I] y_rep[J];
  int M_rep[I];
  real lambda_pred[i];


  // occupancy
  z_pred ~ bernoulli_logit(beta0_psi + beta1_psi * latitude + beta2_psi * longitude + 
                        beta_conifer_psi * percent_conifer + beta_deciduous_psi * percent_deciduous + 
                        beta_pine_psi * percent_pine);
  for (j in 1:J) {
    p_detection_pred[j] ~ beta(beta0_p + beta1_p * time_of_day[j] + beta2_p * Julian_date[j], sigma_ObservationProcess); //
  }

  y_rep ~ Bernoulli(z_pred, p_detection_pred[i]);


  // Abundance
  M_rep ~ poisson_log(z_pred * lambda_pred);  // is this vectorized?
 
  lambda_pred = poisson_log(beta0_lambda + beta_size * size + beta_age * age + 
                            beta_conifer_lambda * percent_conifer + beta_deciduous_lambda * percent_deciduous + 
                            beta_pine_lambda * percent_pine);

} */


generated quantities {
  vector[I] psi_est;  // Estimated occupancy probabilities for each site
  matrix[I, J] p_est;  // Estimated detection probabilities for each site and visit
  vector[I] lambda_est;  // Estimated expected counts for each site when occupied

  for (i in 1:I) {
    // Estimate for occupancy probability for site i
    real lin_pred_occupancy = beta_0 + 
                              beta_latitude * latitude[i] + 
                              beta_longitude * longitude[i] + 
                              beta_percent_conifer * percent_conifer[i] + 
                              beta_percent_deciduous * percent_deciduous[i] +
                              beta_percent_pine * percent_pine[i];
    psi_est[i] = inv_logit(lin_pred_occupancy);
    
    // Estimate for Poisson lambda for site i
    real lin_pred_lambda = gamma_0 + gamma_size * size[i] + gamma_age * age[i];
    lambda_est[i] = exp(lin_pred_lambda);
    
    for (j in 1:J) {
      // Estimate for detection probability at visit j for site i
      real lin_pred_detection = alpha_0 + alpha_1 * time_of_day[i,j] + alpha_2 * Julian_date[i,j];
      p_est[i,j] = exp(lin_pred_detection);
    }
  }
}
