data {
  int N;  // number of samples
  int O;  // origination cohort
  int C;  // classes/group
  int N_unc;
  int N_cen;
  int samp_unc[N_unc];
  int samp_cen[N_cen];
  real<lower=0> dur_unc[N_unc];
  int group_unc[N_unc];
  int cohort_unc[N_unc];
  real occupy_unc[N_unc];
  real env_unc[N_unc];
  real size_unc[N_unc];
  real<lower=0> dur_cen[N_cen];
  int group_cen[N_cen];
  int cohort_cen[N_cen];
  real occupy_cen[N_cen];
  real env_cen[N_cen];
  real size_cen[N_cen];
}
parameters {
  vector[5] mu_prior;
  vector[5] beta[O];  // betas
  corr_matrix[5] Omega;
  vector<lower=0>[5] sigma;
}
transformed parameters {
  cov_matrix[5] Sigma;

  Sigma <- quad_form_diag(Omega, sigma);
}
model {
  // varying slopes, varying intercepts
  // done by cohort
  Omega ~ lkj_corr(2);
  sigma ~ cauchy(0, 1);
  for(i in 1:5) {
    mu_prior[i] ~ normal(0, 5);
  }
  for(i in 1:O) {
    beta[i] ~ multi_normal(mu_prior, Sigma);
  }

  // likelihood
  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      increment_log_prob(exponential_cdf_log(dur_unc[i],
            exp(beta[cohort_unc[i], 1] +
              beta[cohort_unc[i], 2] * occupy_unc[i] +
              beta[cohort_unc[i], 3] * env_unc[i] +
              beta[cohort_unc[i], 4] * (env_unc[i] * env_unc[i]) +
              beta[cohort_unc[i], 5] * size_unc[i])));
    } else {
      increment_log_prob(exponential_log(dur_unc[i],
            exp(beta[cohort_unc[i], 1] +
              beta[cohort_unc[i], 2] * occupy_unc[i] +
              beta[cohort_unc[i], 3] * env_unc[i] +
              beta[cohort_unc[i], 4] * (env_unc[i] * env_unc[i]) +
              beta[cohort_unc[i], 5] * size_unc[i])));
    }
  }
  for(i in 1:N_cen) {
    increment_log_prob(exponential_ccdf_log(dur_cen[i],
          exp(beta[cohort_cen[i], 1] +
            beta[cohort_cen[i], 2] * occupy_cen[i] +
            beta[cohort_cen[i], 3] * env_cen[i] +
            beta[cohort_cen[i], 4] * (env_cen[i] * env_cen[i]) +
            beta[cohort_cen[i], 5] * size_cen[i])));
  }
}
generated quantities {
  vector[N] log_lik;

  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      log_lik[i] <- exponential_cdf_log(dur_unc[i],
          exp(beta[cohort_unc[i], 1] +
            beta[cohort_unc[i], 2] * occupy_unc[i] +
            beta[cohort_unc[i], 3] * env_unc[i] +
            beta[cohort_unc[i], 4] * (env_unc[i] * env_unc[i]) +
            beta[cohort_unc[i], 5] * size_unc[i]));
    } else {
      log_lik[i] <- exponential_log(dur_unc[i],
          exp(beta[cohort_unc[i], 1] +
            beta[cohort_unc[i], 2] * occupy_unc[i] +
            beta[cohort_unc[i], 3] * env_unc[i] +
            beta[cohort_unc[i], 4] * (env_unc[i] * env_unc[i]) +
            beta[cohort_unc[i], 5] * size_unc[i]));
    }
  }
  for(i in 1:N_cen) {
    log_lik[i + N_unc] <- exponential_ccdf_log(dur_cen[i],
        exp(beta[cohort_cen[i], 1] +
          beta[cohort_cen[i], 2] * occupy_cen[i] +
          beta[cohort_cen[i], 3] * env_cen[i] +
          beta[cohort_cen[i], 4] * (env_cen[i] * env_cen[i]) +
          beta[cohort_cen[i], 5] * size_cen[i]));
  }
}
