data {
  int N;
  int T;
  int sight[N, T];
}
parameters {
  real<lower=0,upper=1> p;
  real<lower=0,upper=1> phi;
  real<lower=0,upper=1> gamma;
}
model {
  real run_prob[T];
  real k;
  k <- 1;
  
  run_prob[1] <- gamma;
  for(t in 2:T) {
    k <- k * (1 - run_prob[t - 1]);
    run_prob[t] <- run_prob[t - 1] * phi + 
      gamma * k;
  }
  
  for(n in 1:N) {
    for(t in 1:T) {
      sight[n, t] ~ bernoulli(p * run_prob[t]);
    }
  }
}

