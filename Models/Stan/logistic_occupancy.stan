data {
  int<lower=0> I;
  array[I] int<lower=0, upper=1> y;
    vector[I] age;
}

parameters {

    real alpha;
    real beta_age;
}

model {

  y ~ bernoulli_logit(alpha + beta_age*age); // the stan implementation of the bernoulli_logit, more explicitely, is bernoulli(inv_logit(alpha + beta* x)), where alpha/betas are scalars, x is a vector.

      //priors
    alpha ~ normal(0, 0.5);
    beta_age ~ normal(0,1);
}

generated quantities {
  int y_pred[I];  // predicted counts
  for (i in 1:I) {
    // Generate predicted count for each observation
    y_pred[i] = bernoulli_logit_rng(alpha + beta_age * age[i] );
  }
}