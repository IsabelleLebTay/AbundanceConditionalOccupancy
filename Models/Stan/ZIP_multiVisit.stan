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
    // real<lower=0, upper = 1> psi;

    real alpha_abund; // Intercept
    real beta_age; // Coefficient for the age
    real beta_size;
    real beta_age_size; // Coefficient for the interaction between age and size

    real alpha_theta;
    real beta_latitude;
    real beta_longitude;
}

transformed parameters {
    // dont need this //vector [I] log_lambda;  // Log population size
    matrix[I, J] logit_p;   // Logit detection probability

    for (i in 1:I) {
    // covariates that affect abundance
   // dont want this here // log_lambda[i] = alpha_abund 
                        + beta_age * age[i] 
                        + beta_size * size[i] 
                        + beta_age_size * age[i] * size[i] ;

    logit_p = rep_matrix(alpha_detect 
                        + beta_tod * time_of_day
                        + beta_toy * time_of_year, J);
    for (j in 1:J) {
    // covariates that affect detection probability

      logit_p[i,j] = alpha_detect 
                        + beta_tod* time_of_day[i, j]
                        + beta_toy * time_of_year[i, j]; 
    }
  }
}

model {

    // detectino priors
    alpha_detect ~ normal(0,1);
    beta_tod ~ normal(0,1);
    beta_toy ~ normal(0,1);

    //priors
    alpha_abund ~ normal(0,1);
    beta_age ~ normal(0,1);
    beta_size ~ normal(0,1);
    beta_age_size ~ normal(0,1); // Priors for interaction term

    alpha_theta ~ normal(0,1);
    beta_latitude ~ normal(0, 1);
    beta_longitude ~ normal(0, 1);

    // theta was a parameter:     real<lower=0, upper=1> theta; // probability of occupancy

    // matrix count [I,J] ~ binomial(z[I], psi[I,J])
    // logit(psi[I, J]) = alpha_detect + beta_tod * time_of_day
    // alpha_detect ~ normal(0, 1)


    /* for (i in 1:I){
        for (j in 1:J) {
           real psi = inv_logit(alpha_detect 
                                    + beta_tod * time_of_day[j]
                                    + beta_toy * time_of_year[j]); 

            
            target += binomial_lpmf(count[i, j] | z[i] , psi[i, j]);

            // count[i, j] ~ binomial(z[i], psi[i, j]);
        }
        // this one, theta is not related to a linear equation prior
        //target += binomial_lpmf(count[i, j] | z[i] , psi[i, j]);

        //close:   this is a better marginalising of z.
        real binomial_poisson_lpmf(int[] y, vector p, vector lambda) { // vectorized;  assumes that y, p, and lambda are all of the same length.
    
    real binomial_logit_lpmf(ints n | ints N, reals alpha)
    
    vector[size(y)] out;
    
    for(i in 1:size(y))
      out[i] = -lambda[i] * p[i] + y[i] * (log(p[i]) + log(lambda[i])) - lgamma(y[i]+1));
    return(sum(out));
  }
    } */

    vector [I] theta = inv_logit(alpha_theta + 
                              beta_latitude * latitude + 
                              beta_longitude * longitude);

        // think about what theta means, and where presence/absence = theta  or 1-theta.
        // here, swith to the logit scale


    // feed z into the poisson part as the true unobserved local abundance, instead of M.
    for (i in 1:I) { 

        vector(K - max_count[i] + 1);

        for (j in 1:(K[i] - max_count[i] + 1)) {
                lp[j] = poisson_log_lpmf(max_count[i] + j-1 | log_lambda[i])
                                        + binomial_logit_lpmf(count[i] | max_count[i] + j - 1, logit_p[i]);
                target += log_sum_exp(lp);
        }
        
  

        // zero-inflated Poisson. log likelihood of count, given theta and lambda
        if (z[i] == 0) {
            target += log_sum_exp(log(theta), // log_sum_exp(arg1, arg2) is the same as  log(exp(arg1) + exp(arg2))
                    log1m(theta) // this computes log(1-theta)
                    + poisson_lpmf(0 | log_lambda));
        } 
        
        else {
            target += log1m(theta)
                    + poisson_lpmf(z[i] | log_lambda);
        }
    }
}
