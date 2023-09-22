
/* data {
    int<lower=0> I; 
  int<lower=0> M[I];  // Observed abundance for each site
  
  // site-level covariates for abundance
  vector[I] size;
  vector[I] age;
}

parameters {
  real gamma_0;  // Intercept for abundance
  real gamma_size;
  real gamma_age;
}

model {
  // Priors for abundance coefficients
  gamma_0 ~ normal(0, 5);
  gamma_size ~ normal(0, 5);
  gamma_age ~ normal(0, 5);

  for (i in 1:I) {
    // Linear predictor for Poisson lambda for site i
    real lin_pred_lambda = gamma_0 + gamma_size * size[i] + gamma_age * age[i];
    
    // Poisson likelihood for site i
    M[i] ~ poisson_exp(lin_pred_lambda);
  }
} */

data {
    int<lower=0> I;
    array[I] int<lower=0> M;
}

parameters {
    real<lower=0, upper=1> theta;
    real<lower=0> lambda;
}
model {
    for (i in 1:I) {
        if (M[i] == 0) {
            target += log_sum_exp(log(theta),
                log1m(theta)
                    + poisson_lpmf(M[i] | lambda));
    } else {
        target += log1m(theta)
            + poisson_lpmf(M[i] | lambda);
        }
    }
}