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

    for (i in 1:I) {
        // Calculate lambda using the size covariate
        real lambda_i = exp(alpha + beta_age * age[i] + beta_size * size[i] + beta_age_size * age[i] * size[i] );
        if (M[i] == 0) {
            target += log_sum_exp(log(theta), // log_sum_exp(arg1, arg2) is the same as  log(exp(arg1) + exp(arg2))
                log1m(theta) // this computes log(1-theta)
                    + poisson_lpmf(M[i] | lambda_i));
    } else {
        target += log1m(theta)
            + poisson_lpmf(M[i] | lambda_i);
        }
    }
}

generated quantities {
    real predicted_counts[I];
    for (i in 1:I) {
        predicted_counts[i] = alpha + beta_age * age[i] + beta_size * size[i] + beta_age_size * age[i] * size[i];
    }
}
