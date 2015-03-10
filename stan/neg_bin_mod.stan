data {
  int C;  // total number of occurrences
  int G;  // number of genera
  int O;  // numer of orders
  int count[C];
  int genus[C];
  int order[G];
  real off[C];
}
parameters {
  real<lower=0> mu[G];
  real mu_gen[O];
  real<lower=0> sigma_gen;
  real mu_ord;
  real<lower=0> sigma_ord;

  real<lower=0> phi;
}
transformed parameters {
}
model {

  // classes
  mu_ord ~ normal(0, 10);  // ultimate mean observation rate
  sigma_ord ~ cauchy(0, 2.5);
  for(o in 1:O) {
    mu_gen[o] ~ normal(mu_ord, sigma_ord);
  }
  
  // genera
  sigma_gen ~ cauchy(0, 2.5);
  for(g in 1:G) {
    mu[g] ~ normal(mu_gen[order[g]], sigma_gen);
  }
  // each individual belongs to a genus. classes of genera have different means, 
  // but assumed to have homogeneous variance.

  // individual
  phi ~ cauchy(0, 2.5);
  for(i in 1:C) {
    count[i] ~ neg_binomial_2_log(mu[genus[i]] + log(off[i]), phi);
  }
}
generated quantities {
  vector[C] log_lik;

  for(i in 1:C) {
    log_lik[i] <- neg_binomial_2_log_log(count[i], mu[genus[i]], phi);
  }
}
