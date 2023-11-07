data {
    int<lower=0> I;
    array[I] int<lower=0> M; // this format is newer, int M[I] will be deprecated
    vector[I] age;
    vector[I] size;

    vector[I] percent_conifer;
    vector[I] percent_pine;
    vector[I] percent_deciduous;
    vector[I] neighbour;
}

parameters {
    real beta_0; // Intercept
    real beta_age; // Coefficient for the age
    real beta_size;
//    real beta_age_size; // Coefficient for the interaction between age and size
    real beta_conifer;
    real beta_pine;
    real beta_deciduous;
    real beta_theta;

    real alpha_0;
    real alpha_neighbour;
   // real<lower=0, upper=1> theta; // probability of occupancy
}

model {

    //priors: fix the distributions! THink about the scale of each prior
    beta_0 ~ normal(0,1);
    beta_age ~ normal(0,1);
    beta_size ~ normal(0,1);
//    beta_age_size ~ normal(0,1); // Priors for interaction term
    beta_conifer ~ normal(0, 1);
    beta_pine ~ normal(0, 1);
    beta_deciduous ~ normal(0,1);
    beta_theta ~ normal(0,1);

    alpha_0 ~ normal(0,1);
    alpha_neighbour ~ normal(0,1);
    // theta ~ uniform(0,1);



    // theta was a parameter:     real<lower=0, upper=1> theta; // probability of occupancy



        // think about what theta means, and where prensece/absence = theta  or 1-theta.
        // here, swith to the logit scale
      
    vector [I] theta = inv_logit(alpha_0
                                    + alpha_neighbour * neighbour);

    for (i in 1:I) {

    real lambda = exp(beta_0
                        + beta_age * age[i]
                        + beta_size * size[i]
                        + beta_conifer * percent_conifer[i]
                        + beta_pine * percent_pine[i]
                        + beta_deciduous * percent_deciduous[i]
                        + beta_theta * logit(theta[i]));

        if (M[i] == 0) {
            target += log_sum_exp(log(theta[i]), 
                log1m(theta[i]) // this computes log(1-theta)
                    + poisson_lpmf(M[i] | lambda));
        } 
        
        else {
            target += log1m(theta[i])
                + poisson_lpmf(M[i] | lambda);
        }
    }
}

generated quantities {
    int M_pred[I];
    real log_lik[I]; // log-likelihood for each predicted count
    vector[I] lambda_rep; // to store expected rate parameter
    vector[I] theta_rep;

    for (i in 1:I) {

        theta_rep[i] = inv_logit(alpha_0
                            + alpha_neighbour * neighbour[i]);

            // Compute lambda for current age and size
        lambda_rep[i] = exp(beta_0 
                            + beta_age * age[i] 
                            + beta_size * size[i]
                            + beta_conifer * percent_conifer[i]
                            + beta_pine * percent_pine[i]
                            + beta_deciduous * percent_deciduous[i]
                            + beta_theta * logit(theta_rep[i]));
                            
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
