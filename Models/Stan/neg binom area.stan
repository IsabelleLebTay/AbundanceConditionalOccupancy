data {
    int<lower=0> I; // number of sites
    array[I] int<lower=0> M;
    vector[I] size;
    array[I] int tree_groups;
    vector[I] patch;
    vector[I] age;
    vector[I] age2;
    vector[I] dist_forest;
}

parameters {
    real alpha;
    real beta_size;
    real beta_age_patch;
    real beta_patch;
    vector[4] beta_age;
    real beta_age2;
    real beta_dist_forest;
    real<lower=0> phi; // dispersion parameter for the negative-binomial
}

model {
    vector[I] mu;
    
    for (i in 1:I) {
        mu[i] = exp(alpha 
                    + beta_size * size[i] 
                    + beta_age_patch * patch[i] .* age[i]
                    + beta_age[tree_groups[i]] .* age[i]
                    + beta_patch * patch[i]
                    + beta_dist_forest * dist_forest[i]
                    + beta_age2 * age2[i]);
    }
    
    // likelihood
    M ~ neg_binomial_2(mu, phi);
    
    // priors
    alpha ~ normal(0, 0.5);
    beta_size ~ normal(0, 1);
    beta_patch ~ normal(0, 1);
    beta_age_patch ~ normal(0, 1);
    beta_age ~ normal(0, 1);
    beta_dist_forest ~ normal(0, 1);
    beta_age2 ~ normal(0, 1);
    phi ~ cauchy(0, 5); // weakly informative prior for dispersion parameter
}

generated quantities {
    vector[I] log_lik;
    for (i in 1:I) {
        log_lik[i] = neg_binomial_2_lpmf(M[i] | exp(alpha 
                                                   + beta_size * size[i]
                                                   + beta_age_patch * patch[i] .* age[i]
                                                   + beta_age[tree_groups[i]] .* age[i]
                                                   + beta_patch * patch[i]
                                                   + beta_dist_forest * dist_forest[i]
                                                   + beta_age2 * age2[i]), 
                                         phi);
    }
}
