data {
    int<lower=0> I;
    array[I] int<lower=0> M; // this format is newer, int M[I] will be deprecated
    vector[I] age;
    vector[I] size;

}

parameters {
    real<lower=0, upper=1> theta; // probability of occupancy
    real alpha; // Intercept
    //real beta_age; // Coefficient for the age
    //real beta_size;
    real beta_age_size; // Coefficient for the interaction between age and size
}

model {

    //priors: fix the distributions! THink about the scale of each prior
    alpha ~ normal(0,1);
   // beta_age ~ normal(0,1);
    //beta_size ~ normal(0,1);
    beta_age_size ~ normal(0,1); // Priors for interaction term


    // vector[I] lambda = exp(alpha + beta_age * age + beta_size * size + beta_age_size * age * size );
    // vector[I] lambda = exp(alpha + beta_age * age + beta_size * size);

    for (i in 1:I) {  
        real lambda_i = exp(alpha + beta_age_size * age[i] * size[i] );  
        if (M[i] == 0) {
            target += log_sum_exp(log(theta), // log_sum_exp(arg1, arg2) is the same as  log(exp(arg1) + exp(arg2))
                log1m(theta) // this computes log(1-theta)
                    + poisson_lpmf(M[i] | lambda_i));

        } 
        
        else {
            target += log1m(theta)
                + poisson_lpmf(M[i] | lambda_i);
        }
    }
}
