data {
  int N;  // number of samples
  int O;  // origination cohort
  int T;

  real<lower=0> dur[N];  // duration of uncensored obs

  int censored[N];

  int cohort[N];  // cohort membership
  real occupy[N];  // range size
  real samp[N];  
}
parameters {
  real alpha_trans;

  // regression coefficients
  vector[2] mu_prior;
  vector[2] beta[O];  // cohort coef 
  corr_matrix[2] Omega;
  vector<lower=0>[2] sigma; 

  real delta;
}
transformed parameters {
  real<lower=0> alpha;  // individual shape
  cov_matrix[2] Sigma;

  alpha <- exp(alpha_trans);
  Sigma <- quad_form_diag(Omega, sigma);
}
model {
  real hold[N];

  alpha_trans ~ normal(0, 1);

  // regression coefficients
  Omega ~ lkj_corr(2);
  sigma ~ cauchy(0, 1);
  mu_prior[1] ~ normal(0, 5);
  mu_prior[2] ~ normal(-1, 1);

  delta ~ normal(0, 1);

  for(i in 1:O) {
    beta[i] ~ multi_normal(mu_prior, Sigma);
  }
  
  for(i in 1:N) {
    hold[i] <- exp(-(beta[cohort[i], 1] + 
          beta[cohort[i], 2] * occupy[i] + 
          delta * samp[i]) / alpha);
  }

  // likelihood / sampling statements
  for(i in 1:N) {
    if(censored[i] == 0) {
      increment_log_prob(weibull_log(dur[i], alpha, hold[i]) -
          weibull_ccdf_log(T, alpha, hold[i]));
    } else {
      increment_log_prob(weibull_ccdf_log(dur[i], alpha, hold[i]) -
        weibull_ccdf_log(T, alpha, hold[i]));
    }
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_tilde;
  vector[N] hold;

  for(i in 1:N) {
    hold[i] <- exp(-(beta[cohort[i], 1] +
          beta[cohort[i], 2] * occupy[i] +
          delta * samp[i])/ alpha);
  }

  // log_lik
  for(i in 1:N) {
    if(censored[i] == 0) {
      log_lik[i] <- weibull_log(dur[i], alpha, hold[i]) - 
        weibull_ccdf_log(T, alpha, hold[i]);
    } else {
      log_lik[i] <- weibull_ccdf_log(dur[i], alpha, hold[i]) - 
        weibull_ccdf_log(dur[i], alpha, hold[i]);
    }
  }

  // posterior predictive simulations
  for(i in 1:N) {
    y_tilde[i] <- weibull_rng(alpha, hold[i]);
  }
}

