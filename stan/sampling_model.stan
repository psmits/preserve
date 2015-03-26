data {
  int N;  // number of observations
  int G;  // numer of genera
  int O;  // number of orders
  int count[N];
  int genus[N];
  int order[G];
}
parameters {
  real intercept;
  real mu[G];
  real<lower=0> sigma;

  real nu[O];
  real<lower=0> sigma_a;

  real<lower=0,upper=1> theta;
}
model {
  intercept ~ normal(0, 10);

  for(o in 1:O) {
    nu[o] ~ normal(0, sigma_a);
  }
  sigma_a ~ cauchy(0, 2.5);

  for(g in 1:G) {
    mu[g] ~ normal(nu[order[g]], sigma);
  }
  sigma ~ cauchy(0, 2.5);

  theta ~ beta(1, 1);

  for(i in 1:N) {
    if(count[i] == 0){ 
      increment_log_prob(log_sum_exp(bernoulli_log(1, theta), 
            bernoulli_log(0, theta) + 
            poisson_log_log(count[i], intercept + mu[genus[i]])));
    } else {
      increment_log_prob(bernoulli_log(0, theta) + 
          poisson_log_log(count[i], intercept + mu[genus[i]]));
    }
  }
}
generated quantities {
  vector[C] log_lik;

  for(i in 1:C) {
    log_lik[i] <- poisson_log(count[i], intercept + mu[genus[i]]);
  }
}
