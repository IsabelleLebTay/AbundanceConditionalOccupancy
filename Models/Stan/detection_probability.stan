data{
    int<lower=0> I; // number of sites
    array[I] int<lower=0> y; // detection
    vector[I] size;
    array[I] int tree_groups;
    vector[I] patch;
    vector[I] age;
}

model{
    y ~ 
    p = inv_logit()
}