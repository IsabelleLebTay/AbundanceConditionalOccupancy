data {
    int<lower=0> I;
    array[I] int<lower=0> M; // this format is newer, int M[I] will be deprecated
    vector[I] age;
    vector[I] size;

}

parameters {
    real<lower=0, upper=1> theta; // probability of occupancy
    real alpha; // Intercept
    real beta_age; // Coefficient for the age
    real beta_size;
    real beta_age_size; // Coefficient for the interaction between age and size
}

model {

    //priors: fix the distributions! THink about the scale of each prior
    alpha ~ normal(0,1);
    beta_age ~ normal(0,1);
    beta_size ~ normal(0,1);
    beta_age_size ~ normal(0,1); // Priors for interaction term

//    theta will have to be given a likelihood, work on the logit scale (logistic regression on the extr zeros)
//   think carefully about which side gets theta and which gets 1-theta.


    // linear equation for the abundance parameter. Exp because lambda can only be positive
    // vector[I] lambda = exp(alpha + beta_age * age + beta_size * size + beta_age_size * age * size );
    // vector[I] lambda = exp(alpha + beta_age * age + beta_size * size);
    
    // I haven'ty figured out how to have an interaction term when the linear
    // expression is vectorized this way. It raises on error about matrix operations.
    // Might have to subscript the matrix.

  // for (i in 1:I) {
        // Calculate lambda using the size covariate
    //    real lambda_i = exp(alpha + beta_age * age[i] + beta_size * size[i] + beta_age_size * age[i] * size[i] );
  // }


        // think about what theta means, and where prensece/absence = theta  or 1-theta.
        // here, swith to the logit scale
    for (i in 1:I) {    
        real lambda = exp(alpha 
                            + beta_age * age[i] 
                            + beta_size * size[i] 
                            + beta_age_size * age[i] * size[i] );

        if (M[i] == 0) {
            target += log_sum_exp(log(theta), // log_sum_exp(arg1, arg2) is the same as  log(exp(arg1) + exp(arg2))
                log1m(theta) // this computes log(1-theta)
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


    for (i in 1:I) {

            // Compute lambda for current age and size
        lambda_rep[i] = exp(alpha 
                            + beta_age * age[i] 
                             + beta_size * size[i]
                             + beta_age_size * age[i] * size[i]);
        // Sample from zero-inflated component with probability theta and from Poisson component with probability 1-theta
        if (bernoulli_rng(theta)) {
            M_pred[i] = 0;
            log_lik[i] = log_sum_exp(log(theta), 
                                        log1m(theta) 
                                        + poisson_lpmf(0 | lambda_rep[i]));        
        } 
        
        else {
            M_pred[i] = poisson_rng(lambda_rep[i]);
            log_lik[i] = log1m(theta)
                                + poisson_lpmf(M_pred[i] | lambda_rep[i]);

        }
        
    }
}
