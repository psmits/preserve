data {
  int N;  // number of samples
  int O;  // origination cohort
  int N_unc;  // number uncensored obs
  int N_cen;  // number censored obs
  real<lower=0> dur_unc[N_unc];  // duration of uncensored obs
  int cohort_unc[N_unc];  // cohort membership
  real occupy_unc[N_unc];  // range size
  real size_unc[N_unc];  // body size
  real env_unc[N_unc];  // number of epicontinental occ
  real samp_unc[N_unc];  // number of open ocean occ
  real<lower=0> dur_cen[N_cen];
  int cohort_cen[N_cen];
  real occupy_cen[N_cen];
  real env_cen[N_cen];
  real samp_cen[N_cen];
  real size_cen[N_cen];
}
parameters {
  real alpha_trans;

  // regression coefficients
  vector[8] mu_prior;
  vector[8] beta[O];  // cohort coef 
  corr_matrix[8] Omega;
  vector<lower=0>[8] sigma; 
}
transformed parameters {
  real<lower=0> alpha;  // individual shape
  cov_matrix[8] Sigma;

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
  mu_prior[6] ~ normal(-1, 1);
  mu_prior[7] ~ normal(0, 1);
  mu_prior[8] ~ normal(0, 1);

  for(i in 1:O) {
    beta[i] ~ multi_normal(mu_prior, Sigma);
  }

  // likelihood / sampling statements
  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      increment_log_prob(weibull_cdf_log(dur_unc[i], alpha,
            exp(-(beta[cohort_unc[i], 1] + 
                beta[cohort_unc[i], 2] * occupy_unc[i] + 
                beta[cohort_unc[i], 3] * env_unc[i] + 
                beta[cohort_unc[i], 4] * (env_unc[i]^2) +
                beta[cohort_unc[i], 5] * size_unc[i] +
                beta[cohort_unc[i], 6] * samp_unc[i] +
                beta[cohort_unc[i], 7] * (occupy_unc[i] * samp_unc[i]) +
                beta[cohort_unc[i], 8] * (env_unc[i] * samp_unc[i])) 
              / alpha)));
    } else {
      increment_log_prob(weibull_log(dur_unc[i], alpha,
            exp(-(beta[cohort_unc[i], 1] + 
                beta[cohort_unc[i], 2] * occupy_unc[i] + 
                beta[cohort_unc[i], 3] * env_unc[i] + 
                beta[cohort_unc[i], 4] * (env_unc[i]^2) +
                beta[cohort_unc[i], 5] * size_unc[i] +
                beta[cohort_unc[i], 6] * samp_unc[i] +
                beta[cohort_unc[i], 7] * (occupy_unc[i] * samp_unc[i]) +
                beta[cohort_unc[i], 8] * (env_unc[i] * samp_unc[i])) 
              / alpha)));
    }
  }
  for(i in 1:N_cen) {
    increment_log_prob(weibull_ccdf_log(dur_unc[i], alpha,
          exp(-(beta[cohort_cen[i], 1] + 
              beta[cohort_cen[i], 2] * occupy_cen[i] + 
              beta[cohort_cen[i], 3] * env_cen[i] + 
              beta[cohort_cen[i], 4] * (env_cen[i]^2) +
              beta[cohort_cen[i], 5] * size_cen[i] +
              beta[cohort_cen[i], 6] * samp_cen[i] +
              beta[cohort_cen[i], 7] * (occupy_cen[i] * samp_cen[i]) +
              beta[cohort_cen[i], 8] * (env_cen[i] * samp_cen[i])) 
            / alpha)));
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_tilde;
  vector[N] hold;

  for(i in 1:N_unc) {
    hold[i] <- exp(-(beta[cohort_unc[i], 1] +
          beta[cohort_unc[i], 2] * occupy_unc[i] +
          beta[cohort_unc[i], 3] * env_unc[i] + 
          beta[cohort_unc[i], 4] * (env_unc[i]^2) +
          beta[cohort_unc[i], 5] * size_unc[i] +
          beta[cohort_unc[i], 6] * samp_unc[i] +
          beta[cohort_unc[i], 7] * (occupy_unc[i] * samp_unc[i]) +
          beta[cohort_unc[i], 8] * (env_unc[i] * samp_unc[i])) / alpha);
  }
  for(i in 1:N_cen) {
    hold[i + N_unc] <- exp(-(beta[cohort_cen[i], 1] +
          beta[cohort_cen[i], 2] * occupy_cen[i] +
          beta[cohort_cen[i], 3] * env_cen[i] + 
          beta[cohort_cen[i], 4] * (env_cen[i]^2) +
          beta[cohort_cen[i], 5] * size_cen[i] +
          beta[cohort_cen[i], 6] * samp_cen[i] +
          beta[cohort_cen[i], 7] * (occupy_cen[i] * samp_cen[i]) +
          beta[cohort_cen[i], 8] * (env_cen[i] * samp_cen[i])) / alpha);
  }

  // log_lik
  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      log_lik[i] <- weibull_cdf_log(dur_unc[i], alpha, hold[i]);
    } else {
      log_lik[i] <- weibull_log(dur_unc[i], alpha, hold[i]);
    }
  }
  for(i in 1:N_cen) {
    log_lik[i + N_unc] <- weibull_ccdf_log(dur_cen[i], 
        alpha, hold[i]);
  }

  // posterior predictive simulations
  for(i in 1:N_unc) {
    y_tilde[i] <- weibull_rng(alpha, hold[i]);
  }
  for(i in 1:N_cen) {
    y_tilde[i + N_unc] <- weibull_rng(alpha, hold[i]);
  }
}

