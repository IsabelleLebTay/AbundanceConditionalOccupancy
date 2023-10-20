data {
    int<lower=0> I;            // Number of sites
    array[I] int<lower=0> M; // this format is newer, int M[I] will be deprecated
    int<lower=0> J;              // Number of visits
    matrix[I,J] count;            // Observed counts per site and visit
    matrix[I, J] time_of_day;
    matrix[I, J] time_of_year;

    vector[I] age;
    vector[I] size;

    vector[I] latitude;
    vector[I] longitude;
}

parameters {
    real alpha_detect;
    real beta_tod;
    real beta_toy;

    real alpha_abund; // Intercept
    real beta_age; // Coefficient for the age
    real beta_size;
    real beta_age_size; // Coefficient for the interaction between age and size

    real alpha_theta;
    real beta_latitude;
    real beta_longitude;
}

model {

    // detectino priors
    alpha_detect ~ normal(0,1)

    //priors
    alpha_abund ~ normal(0,1);
    beta_age ~ normal(0,1);
    beta_age2 ~ normal(0,1);
    beta_size ~ normal(0,1);
    beta_age_size ~ normal(0,1); // Priors for interaction term

    alpha_theta ~ normal(0,1);
    beta_latitude ~ normal(0, 1);
    beta_longitude ~ normal(0, 1);

    // theta was a parameter:     real<lower=0, upper=1> theta; // probability of occupancy

    // matrix count [I,J] ~ binomial(z[I], psi[I,J])
    // logit(psi[I, J]) = alpha_detect + beta_tod * time_of_day
    // alpha_detect ~ normal(0, 1)


    for (i in 1:I){
        for (j in 1:J) {
            logit(psi[i, j]) = alpha_detect 
                                + beta_tod * time_of_day[j]
                                + beta_toy * time_of_year[j];

            count[i, j] ~ binomial(z[i], psi[i, j])

        }
    }

    vector [I] theta = inv_logit(alpha_theta + 
                              beta_latitude * latitude + 
                              beta_longitude * longitude);

        // think about what theta means, and where presence/absence = theta  or 1-theta.
        // here, swith to the logit scale


    // feed z into the poisson part as the true unobserved local abundance, instead of M.
    for (i in 1:I) {    
        real lambda = exp(alpha_abund 
                            + beta_age * age[i] 
                            + beta_size * size[i] 
                            + beta_age_size * age[i] * size[i] );

        if (z[i] == 0) {
            target += log_sum_exp(log(theta[i]), // log_sum_exp(arg1, arg2) is the same as  log(exp(arg1) + exp(arg2))
                log1m(theta[i]) // this computes log(1-theta)
                    + poisson_lpmf(z[i] | lambda));

        } 
        
        else {
            target += log1m(theta)
                + poisson_lpmf(z[i] | lambda);
        }
    }
}
