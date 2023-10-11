// this is the Stan ZIP implementation in the documentation

data {
    int<lower=0> I;
    array[I] int<lower=0> M;
    vector[I] age;
    vector[I] size;
//    vector[I] percent_conifer;
}


// we need to both tell the Stan compiler an undefined function is okay and let C++ know what it should be.
functions {
    int num_zeros(array[] int M) {
        int sum = 0;
        for (n in 1:size(M)) {
            sum += (M[n] == 0);
            }
    return sum;
    }
}


// this block explicitely separestes the instnace of 0 and non-zeros into two groups
transformed data {
    int<lower=0> N_zero = num_zeros(M);
    array[N - N_zero] int<lower=1> M_nonzero;
    int N_nonzero = 0;
    for (n in 1:N) {
        if (M[n] == 0) continue;
        N_nonzero += 1;
        M_nonzero[N_nonzero] = M[n];
    }
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
    beta_age ~ normal(0,1);
    beta_size ~ normal(0,1);
//    beta_conifer ~ normal(0,1);
    beta_age_size ~ normal(0,1); // Priors for interaction term

    // not vectorized
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

model {
    // zeros first:
    target
        += N_zero
            * log_sum_exp(log(theta),
                            log1m(theta)
                            + poisson_lpmf(0 | lambda));

    // counts next:
    target += N_nonzero * log1m(theta);
    target += poisson_lpmf(y_nonzero | lambda);
   
}

/*
generated quantities {
    real predicted_counts[I];
    for (i in 1:I) {
        predicted_counts[i] = alpha + beta_age * age[i] + beta_size * size[i] + beta_age_size * age[i] * size[i];
    }
}
*/
