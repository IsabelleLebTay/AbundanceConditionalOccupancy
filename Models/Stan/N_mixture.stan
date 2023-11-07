data {
        int<lower=0> I;                 // Number of sites
        int<lower=0> J;                 // Number of visits
        array[I] int<lower=0> count;    // count per site&visit
        
        matrix[R, T] temp;
        int<lower=0> K;                 // upper bound of popualtion size
}

transformed data {
    int<lower=0> max_count[I];

    for (i in 1:I) {
        max_count[i] = max(count[i]);
    }
}

parameter {
    real alpha0;
    real alpha1;
    rela beta0;
    real beta1;
}

transformed parameter {
    vector[I] log_lambda;   // Log population size
    matrix[I, J] logit_p;   // Logit detection probability

    log_lambda = alpha0 + alpha1 * age;
    logit_p = rep_matrix(beta0 + beta1 * age, J);
}

model {
    // Priors

    // Likelihood
    for (i in 1:I) {
        vector(K - max_count[i] + 1)

        for (j in 1(K - max_count[i] +1))
        lp[j] = poisson_log_lpmf(max_count[i] + j-1 | log_lambda[i])
                + binomial_logit_lpmf(count[i] | max_count[i] + j - 1, logit_p[i]);
        target += log_sum_exp(lp);
    }
}

generated quantities {
    int N[I];
    int totalN;

    for (i in 1:I) {
        N[i] = poisson_log_rng(log_lambda[i]);
    totalN = sum(N)
    }
}