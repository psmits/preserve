data {
  int N;  // number of samples
  int O;  // origination cohort
  int K;  // number covariates
  real<lower=0> y[N];  // duration
  matrix[N, K] X;

  int status[N];  // 1 == censored, 0 == uncensored
  int cohort[N];
}
parameters {
  real alpha_trans;

  vector[K] beta;

  vector[O] c;
  real<lower=0> sigma_c;
}
transformed parameters {
  real<lower=0> alpha;
  vector<lower=0>[N] sigma;

  // shape parameter
  alpha = exp(alpha_trans);
  
  // scale parameter
  sigma[1:N] = exp(X[1:N] * beta + c[cohort]);
}
model {
  alpha_trans ~ normal(0, 1);

  beta ~ normal(0, 5);

  c ~ normal(0, sigma_c);
  sigma_c ~ normal(0, 5);

  for(i in 1:N) {
    if(status[i] == 0) {
      target += weibull_lpdf(y[i] | alpha, sigma[i]);
    } else {
      target += weibull_lccdf(y[i] | alpha, sigma[i]);
    }
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_tilde;

  // log_lik
  for(i in 1:N) {
    if(status[i] == 0) {
      log_lik[i] = weibull_lpdf(y[i] | alpha, sigma[i]);
    } else {
      log_lik[i] = weibull_lccdf(y[i] | alpha, sigma[i]);
    }
  }

  // posterior predictive simulations
  for(i in 1:N) {
    y_tilde[i] = weibull_rng(alpha, sigma[i]);
  }
}
