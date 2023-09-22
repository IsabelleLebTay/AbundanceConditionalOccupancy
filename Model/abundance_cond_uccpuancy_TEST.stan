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
  int<lower=0> y[I, J];  // observed presence/absence, derived from the count N
  real<lower=0> M[I];       // mean bird count for each site

  for (i in 1:I) {

    M[i] = sum(N[i,]) / J;  // compute the mean count for site i

    for (j in 1:J) {
      y[i, j] = (N[i, j] > 0);
    }
  }
}


parameters {
  // Occupancy
  vector[I]<lower=0> p_detection[J]; // the probability of detection is vecotrized over the number of sites for length of number of visits
  int<lower=0> z[I]; // int z, estimated occupancy, of size I (0 or 1)
  real beta0_psi;  // intercept for occupancy
  real beta1_psi;  // coefficient for latitude
  real beta2_psi;  // coefficient for longitude
  real beta_conifer_psi;
  real beta_deciduous_psi;
  real beta_pine_psi;

  // Detection
  real beta0_p;    // intercept for detection probability
  real beta1_p;    // coefficient for time of day
  real beta2_p;    // coefficient for Julian date
  real sigma_ObservationProcess;

  // Abundance
  real beta0_lambda;     // intercept for abundance
  real beta_size;  // coefficient for size
  real beta_age;   // coefficient for age
  real beta_conifer_lambda;
  real beta_deciduous_lambda;
  real beta_pine_lambda;
}


model {
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
}

generated quantities{
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

}