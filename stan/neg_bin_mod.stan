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
  mu_ord ~ normal(0, 10);
  sigma_ord ~ cauchy(0, 2.5);
  for(o in 1:O) {
    mu_gen[o] ~ normal(mu_ord, sigma_ord);
  }
  
  sigma_gen ~ cauchy(0, 2.5);
  for(g in 1:G) {
    mu[g] ~ normal(mu_gen[order[g]], sigma_gen);
  }

  phi ~ cauchy(0, 2.5);
  for(i in 1:C) {
    count[i] ~ neg_binomial_2_log(mu[genus[i]] + log(off[i]), phi);
  }
}
