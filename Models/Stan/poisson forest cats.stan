data{
    int<lower=0> I; // number of sites
    array[I] int<lower=0> M;
    vector[I] size;
    array[I] int tree_groups;
}

parameters{
    real alpha;
    vector[4] beta_size;
}

model{

    // likelihood
    M ~ poisson_log(alpha + beta_size[tree_groups] .* size);
    //priors
    alpha ~ normal(0, 0.5);
    beta_size ~ normal(0,1);
}