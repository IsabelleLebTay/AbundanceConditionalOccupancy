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
    vector[4] beta_size;
    real beta_age_patch;

}

model{

    // likelihood
    M ~ poisson_log(alpha + beta_size[tree_groups] .* size
                          + beta_age_patch * patch .* age);

    //priors
    alpha ~ normal(0, 0.5);
    beta_size ~ normal(0,1);
    beta_age_patch ~ normal(0,1);
} 

/*
generated quantities{
    // # how does lambda change in the 4 different habitat types, as the patch size increases?
    // # how does lambda change over age of harvest when patch is either 0 or 1?
    //vector[4] pred_lambda[I];
    array[I] real<lower=0> pred_lambda;
    //real<lower=0> pred_lambda;

    // Generate lambda for different habitat types and sizes

    for(i in 1:I){

        pred_lambda = exponential_rng(alpha + beta_size[tree_groups] .* size[i]
                            + beta_age_patch * patch[i] .* age[i]);
    }

}*/

