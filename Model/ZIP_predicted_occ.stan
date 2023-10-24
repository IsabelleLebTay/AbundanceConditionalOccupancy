data {
    int<lower=0> I;
    array[I] int<lower=0> M; // this format is newer, int M[I] will be deprecated
    vector[I] age;
    vector[I] size;
    vector[I] latitude;
    vector[I] longitude;
    vector[I] percent_conifer;
    vector[I] percent_pine;


}

parameters {
    real alpha; // Intercept
    real beta_age; // Coefficient for the age
    real beta_size;
    real beta_age_size; // Coefficient for the interaction between age and size

    real alpha_theta;
    real beta_latitude;
    real beta_longitude;

    real beta_conifer;
    real beta_pine;
}

model {

    //priors: fix the distributions! THink about the scale of each prior
    alpha ~ normal(0,1);
    beta_age ~ normal(0,1);
    beta_size ~ normal(0,1);
    beta_age_size ~ normal(0,1); // Priors for interaction term

    alpha_theta ~ normal(0,1);
    beta_latitude ~ normal(0, 1);
    beta_longitude ~ normal(0, 1);

    beta_conifer ~ normal(0, 1);
    beta_pine ~ normal(0, 1);


    // theta was a parameter:     real<lower=0, upper=1> theta; // probability of occupancy

    vector [I] theta = inv_logit(alpha_theta + 
                              beta_latitude * latitude + 
                              beta_longitude * longitude);

        // think about what theta means, and where prensece/absence = theta  or 1-theta.
        // here, swith to the logit scale
    for (i in 1:I) {    
        real lambda = exp(alpha 
                            + beta_age * age[i] 
                            + beta_size * size[i] 
                            + beta_age_size * age[i] * size[i] 
                            + beta_conifer * percent_conifer[i]
                            + beta_pine * percent_pine[i]);

        if (M[i] == 0) {
            target += log_sum_exp(log(theta[i]), // log_sum_exp(arg1, arg2) is the same as  log(exp(arg1) + exp(arg2))
                log1m(theta[i]) // this computes log(1-theta)
                    + poisson_lpmf(M[i] | lambda));

        } 
        
        else {
            target += log1m(theta)
                + poisson_lpmf(M[i] | lambda);
        }
    }
}

generated quantities {
    int M_pred[I]; // change this to M_pred
    real log_lik[I]; // log-likelihood for each predicted count
    vector[I] lambda_rep; // to store expected rate parameter
    vector[I] theta_rep;

    for (i in 1:I) {
        
        theta_rep[i] = inv_logit(alpha_theta + 
                              beta_latitude * latitude[i] + 
                              beta_longitude * longitude[i]);

            // Compute lambda for current age and size
        lambda_rep[i] = exp(alpha 
                            + beta_age * age[i] 
                            + beta_size * size[i]
                            + beta_age_size * age[i] * size[i]
                            + beta_conifer * percent_conifer[i]
                            + beta_pine * percent_pine[i]);
        // Sample from zero-inflated component with probability theta and from Poisson component with probability 1-theta
        if (bernoulli_rng(theta_rep[i])) {
            M_pred[i] = 0;
            log_lik[i] = log_sum_exp(log(theta_rep[i]), 
                                        log1m(theta_rep[i]) 
                                        + poisson_lpmf(0 | lambda_rep[i]));        
        } 
        
        else {
            M_pred[i] = poisson_rng(lambda_rep[i]);
            log_lik[i] = log1m(theta_rep[i])
                                + poisson_lpmf(M_pred[i] | lambda_rep[i]);

        }
        
    }
}
