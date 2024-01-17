data {
  int<lower=0> I; //Number of Sites
    int<lower=6> J; //Number of Replicates at each site
    int<lower=0, upper=1> y[I, J]; //Detection at each site on each sampling rep
    int<lower=0, upper=1> x[I]; //Observed occupancy at each site
    real time_of_day[I, J];
    real Julian_date[I, J];
    }
    
parameters {
    real a0; //specifying regression parameters
    real b0;
    real b1;
    real b2;
    }
    
transformed parameters {
    real<lower=0,upper=1>  psi[I]; 
    real<lower=0,upper=1>  p[I, J];
    for(i in 1:I) {
    psi[i] = inv_logit(a0);  //intercept-only model for occupancy
    for(j in 1:J) {
    p[i, j] = inv_logit(b0 + b1*time_of_day[i, j] + b2*Julian_date[i, j]);
    }
    }
}
    
model {
    // Priors
    a0 ~ normal(0, 1);
    b0 ~ normal(0, 1);
    b1 ~ normal(0, 1);
    b2 ~ normal(0, 1);
    
    // likelihood
    for(i in 1:I) {
        if(x[i]==1) {
            1 ~ bernoulli(psi[i]);
            y[i] ~ bernoulli(p[i]);
            }
        if(x[i] == 0) {
            real log_prob_not_detected = 0;
            for(j in 1:J) {
                log_prob_not_detected += log1m(p[i, j]);
            }
            increment_log_prob(log_sum_exp(log(psi[i]) + log_prob_not_detected, log1m(psi[i])));
        }
    }
}
