data{
    int<lower=0> I; // number of sites
    array[I] int<lower=0> M;
    vector[I] edge;
    array[I] int tree_groups;
    vector[I] patch;
    vector[I] age;
    vector[I] dist_forest;
}

parameters{
    real alpha;
    //vector[4] beta_size; use this for the original formula
    real beta_edge;
    //real beta_age_size;
    real beta_age_patch;
    real beta_patch;
    vector[4] beta_age;
    real beta_dist_forest;
}

model{

    // likelihood


    // in this one, oatch sie doesn't interact with forest type, because ecologically it doesnt make much sense.
    M ~ poisson_log(alpha 
                        + beta_edge * edge
                        // + beta_age_size * size .* age // .* means element-wise multiplication
                        + beta_age_patch * patch .* age
                        + beta_age[tree_groups] .* age
                        + beta_patch * patch
                        + beta_dist_forest * dist_forest); //
    //priors
    alpha ~ normal(0, 0.5);
    beta_edge ~ normal(0,1);
    //beta_age_size ~ normal(0,1);
    beta_patch ~ normal(0,1);
    beta_age_patch ~ normal(0,1);
    beta_age ~ normal(0,1);
    beta_dist_forest ~ normal(0,1);
} 

generated quantities {
  vector[I] log_lik;  
  //real log_lik[I]; // log-likelihood for each predicted count
  for (i in 1:I) {
    // Generate likelihood for each observation

    log_lik[i] = poisson_log_lpmf(M[i] | alpha 
                                        + beta_edge * edge[i]
                                        + beta_age_patch * patch[i] .* age[i]
                                        + beta_age[tree_groups[i]] .* age[i]
                                        + beta_patch * patch[i]
                                        + beta_dist_forest * dist_forest[i]);
  }
}
