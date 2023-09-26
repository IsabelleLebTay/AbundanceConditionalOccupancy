
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
    vector[I] age;
    vector[I] size;
//    vector[I] percent_conifer;
}

parameters {
    real<lower=0, upper=1> theta;
    real alpha; // Intercept
    real beta_age; // Coefficient for the age
    real beta_size;
//    real beta_conifer;
    real beta_age_size; // Coefficient for the interaction between age and size
}

model {
    beta_age ~ normal(0,10);
    beta_size ~ normal(0,10);
//    beta_conifer ~ normal(0,10);
    beta_age_size ~ normal(0,10); // Priors for interaction term

    theta will have to be given a likelihood, work on the logit scale (logistic regression on the extr zeros)
    think carefully about which side gets theta and which gets 1-theta.

    

vector[I] lambda = exp(alpha + beta_age * age + beta_size * size + beta_age_size * age * size );

    for (i in 1:I) {
        // Calculate lambda using the size covariate
        real lambda_i = exp(alpha + beta_age * age[i] + beta_size * size[i] + beta_age_size * age[i] * size[i] );
        if (M[i] == 0) {
            target += log_sum_exp(log(theta), // log_sum_exp(arg1, arg2) is the same as  log(exp(arg1) + exp(arg2))
                log1m(theta) // this computes log(1-theta)
                    + poisson_lpmf(M[i] | lambda[I]));

    } else {
        target += log1m(theta)
            + poisson_lpmf(M[i] | lambda_i);
        }
    }
}

generated quantities {
    real predicted_counts[I];

    you don't have to give priors or declare lambda_pred or theta_pred

    have to recalcuta lambda[i]
    then simulate presence/absence, for all the 1,s simualte abundance

vector[I] lambda = exp(alpha + beta_age * age + beta_size * size + beta_age_size * age * size );
for loop that simulates 0 and 1 as a vernoulli with theta 
there is the proba that bird is absent

    for (i in 1:I) {
        predicted_counts[i] = alpha + beta_age * age[i] + beta_size * size[i] + beta_age_size * age[i] * size[i];
    }
}