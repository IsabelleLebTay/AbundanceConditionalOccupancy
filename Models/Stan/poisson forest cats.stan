data{
    int<lower=0> I; // number of sites
    array[I] int<lower=0> M;
    vector[I] size;
    array[I] int tree_groups;
    vector[I] patch;
    vector[I] age;
}

parameters{
    real alpha;
    //vector[4] beta_size; use this for the original formula
    real beta_size;
    //real beta_age_size;
    real beta_age_patch;
    real beta_patch;
    vector[4] beta_age;
}

model{

    // likelihood
    /* M ~ poisson_log(alpha + beta_size[tree_groups] .* size
                          // + beta_age_size * size .* age // .* means element-wise multiplication
                          + beta_age_patch * patch .* age
                          + beta_age[tree_groups] .* age
                          + beta_patch * patch); */ //this is the original

    M ~ poisson_log(alpha + beta_size * size
                        // + beta_age_size * size .* age // .* means element-wise multiplication
                        + beta_age_patch * patch .* age
                        + beta_age[tree_groups] .* age
                        + beta_patch * patch); //just checking if there is less uncertainty with size when its not per forest type
    //priors
    alpha ~ normal(0, 0.5);
    beta_size ~ normal(0,1);
    //beta_age_size ~ normal(0,1);
    beta_patch ~ normal(0,1);
    beta_age_patch ~ normal(0,1);
    beta_age ~ normal(0,1);
} 


