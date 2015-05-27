data {
  int N;
  int T;
  int sight[N, T];
}
parameters {
  vector[T] p_norm;
  vector[T - 1] phi_norm;
  vector[T] gamma_norm;
  real<lower=0> scales[3];
}
transformed parameters {
  vector<lower=0,upper=1>[T] p;
  vector<lower=0,upper=1>[T - 1] phi;
  vector<lower=0,upper=1>[T] gamma;

  for(t in 1:T) {
    p[t] <- inv_logit(p_norm[t]);
    if(t < T) {
      phi[t] <- inv_logit(phi_norm[t]);
    }
    gamma[t] <- inv_logit(gamma_norm[t]);
  }
}
model {
  real run_prob[T];
  real k;
  k <- 1;
  
  run_prob[1] <- gamma[1];
  for(t in 2:T) {
    k <- k * (1 - run_prob[t - 1]);
    run_prob[t] <- run_prob[t - 1] * phi[t - 1] + 
      gamma[t] * k;
  }
  
  scales[1] ~ cauchy(0, 1);
  scales[2] ~ cauchy(0, 1);
  scales[3] ~ cauchy(0, 1);
  for(t in 1:T) {
    p_norm[t] ~ normal(0, scales[1]);
    if(t < T) {
      phi_norm[t] ~ normal(0, scales[2]);
    }
    gamma_norm[t] ~ normal(0, scales[3]);
  }
  
  for(n in 1:N) {
    for(t in 1:T) {
      sight[n, t] ~ bernoulli(p[t] * run_prob[t]);
    }
  }
}
