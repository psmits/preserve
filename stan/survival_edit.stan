data {
  int N;  // number of samples
  int O;  // origination cohort

  real<lower=0> dur[N];  // duration of uncensored obs

  int censored[N];
  int inclusion[N];

  int cohort[N];  // cohort membership
  real occupy[N];  // range leng
  real env[N];  // number of epicontinental occ
  real leng[N];  // body leng
}
parameters {
  real alpha_trans;

  // regression coefficients
  vector[5] mu_prior;
  vector[5] beta[O];  // cohort coef 
  corr_matrix[5] Omega;
  vector<lower=0>[5] sigma; 
}
transformed parameters {
  real<lower=0> alpha;  // individual shape
  cov_matrix[5] Sigma;

  alpha <- exp(alpha_trans);
  Sigma <- quad_form_diag(Omega, sigma);
}
model {
  alpha_trans ~ normal(0, 1);

  // regression coefficients
  Omega ~ lkj_corr(2);
  sigma ~ cauchy(0, 1);
  mu_prior[1] ~ normal(0, 5);
  mu_prior[2] ~ normal(-1, 1);
  mu_prior[3] ~ normal(0, 1);
  mu_prior[4] ~ normal(1, 1);
  mu_prior[5] ~ normal(0, 1);


  for(i in 1:O) {
    beta[i] ~ multi_normal(mu_prior, Sigma);
  }

  // likelihood / sampling statements
  for(i in 1:N) {
    if(censored[i] == 0) {
      if(dur[i] == 1) {
        increment_log_prob(weibull_cdf_log(dur[i], alpha,
              exp(-(beta[cohort[i], 1] + 
                  beta[cohort[i], 2] * occupy[i] + 
                  beta[cohort[i], 3] * env[i] + 
                  beta[cohort[i], 4] * (env[i]^2) +
                  beta[cohort[i], 5] * leng[i]) / alpha)));
      } else {
        increment_log_prob(weibull_log(dur[i], alpha,
              exp(-(beta[cohort[i], 1] + 
                  beta[cohort[i], 2] * occupy[i] + 
                  beta[cohort[i], 3] * env[i] + 
                  beta[cohort[i], 4] * (env[i]^2) +
                  beta[cohort[i], 5] * leng[i]) / alpha)));
      }
    } else {
      increment_log_prob(weibull_ccdf_log(dur[i], alpha,
            exp(-(beta[cohort[i], 1] + 
                beta[cohort[i], 2] * occupy[i] + 
                beta[cohort[i], 3] * env[i] + 
                beta[cohort[i], 4] * (env[i]^2) +
                beta[cohort[i], 5] * leng[i]) / alpha)));
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
          beta[cohort[i], 3] * env[i] + 
          beta[cohort[i], 4] * (env[i]^2) +
          beta[cohort[i], 5] * leng[i]) / alpha);
  }

  // log_lik
  for(i in 1:N) {
    if(censored[i] == 0) {
      if(dur[i] == 1) {
        log_lik[i] <- weibull_cdf_log(dur[i], alpha, hold[i]);
      } else {
        log_lik[i] <- weibull_log(dur[i], alpha, hold[i]);
      }
    } else {
      log_lik[i] <- weibull_ccdf_log(dur[i], 
          alpha, hold[i]);
    }
  }

  // posterior predictive simulations
  for(i in 1:N) {
    y_tilde[i] <- weibull_rng(alpha, hold[i]);
  }
}
