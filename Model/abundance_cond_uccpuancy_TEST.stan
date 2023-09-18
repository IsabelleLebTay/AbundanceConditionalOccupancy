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
  real M[I];       // mean bird count for each site

  for (i in 1:I) {

    M[i] = sum(N[i,]) / J;  // compute the mean count for site i

    for (j in 1:J) {
      y[i, j] = (N[i, j] > 0);
    }
  }
}


parameters {
  vector[I] p_detection[J]; // the probability of detection is vecotrized over the number of sites for length of number of visits
  real z[I]; // float z, estimated occupancy, of size I 
  real 




}


model {
  //maths
  // z ~ bernoulli_logit(psi, sigma_occupancy_process); //  estimated latent occupancy
  // y ~ Bernoulli(z , p_detection); // vectorized, the likelihood of the observed data of presence/absence

  // not vectorized
  //for i in (1:I) {
  //  z[i] ~ bernoulli_logit(beta0_psi + beta1_psi * latitude[i] + beta2_psi * longitude[i] + 
	//													beta_conifer_psi * percent_conifer[i] + beta_deciduous_psi * percent_deciduous[i] + 
	//													beta_pine_psi * percent_pine[i], sigma_occupancy_process[i]);

  //  for j in (1:J) {
  //    p_detection[i, j] ~ bernoulli_logit(beta0_p + beta1_p * time_of_day[i, j] + beta2_p * Julian_date[i, j]); //
  //  }

  // vectorized likelihoods
  // occupancy
    z ~ bernoulli_logit(beta0_psi + beta1_psi * latitude + beta2_psi * longitude + 
                          beta_conifer_psi * percent_conifer + beta_deciduous_psi * percent_deciduous + 
                          beta_pine_psi * percent_pine, sigma_occupancy_process);
      for j in (1:J) {
        p_detection[j] ~ beta0_p + beta1_p * time_of_day[j] + beta2_p * Julian_date[j] //
      }
  // abundance, not vectorized
//    for i in (1:I) {
//      M[i] ~ poisson_log(z[i] * (beta0_lambda + beta_size * size[i] + beta_age * age[i] + 
//														 beta_conifer_lambda * percent_conifer[i] + beta_deciduous_lambda * percent_deciduous[i] + 
//														 beta_pine_lambda * percent_pine[i]))
//    }

    M ~ poisson_log(z * lambda);
    lambda ~


    // Priors
    sigma_occupancy_process ~ cauchy(0,5);
    

  }
   


  alpha ~ normarl(0, 10);
  beta ~ normal(0, 10);
  sigma ~ cauchy(0,5);

}

generated quantities{



}