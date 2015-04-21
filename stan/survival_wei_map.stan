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
  real lit_unc[N_unc];
  real size_unc[N_unc];
  real<lower=0> dur_cen[N_cen];
  int group_cen[N_cen];
  int cohort_cen[N_cen];
  real occupy_cen[N_cen];
  real env_cen[N_cen];
  real lit_cen[N_cen];
  real size_cen[N_cen];
}
parameters {
  real<lower=0> alpha;
  vector[6] mu_prior;
  vector[6] beta[O];  // betas
  corr_matrix[6] Omega;
  vector<lower=0>[6] sigma;
}
transformed parameters {
  cov_matrix[6] Sigma;

  Sigma <- quad_form_diag(Omega, sigma);
}
model {

  alpha ~ cauchy(0, 2);
  // varying slopes, varying intercepts
  // done by cohort
  Omega ~ lkj_corr(2);
  sigma ~ cauchy(0, 1);
  for(i in 1:6) {
    mu_prior[i] ~ normal(0, 5);
  }
  for(i in 1:O) {
    beta[i] ~ multi_normal(mu_prior, Sigma);
  }

  // likelihood
  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      increment_log_prob(weibull_cdf_log(dur_unc[i], alpha, 
            exp(-(beta[cohort_unc[i], 1] +
                beta[cohort_unc[i], 2] * occupy_unc[i] +
                beta[cohort_unc[i], 3] * env_unc[i] +
                beta[cohort_unc[i], 4] * lit_unc[i] +
                beta[cohort_unc[i], 5] * (env_unc[i] * lit_unc[i]) +
                beta[cohort_unc[i], 6] * size_unc[i]) / alpha)));
    } else {
      increment_log_prob(weibull_log(dur_unc[i], alpha,
            exp(-(beta[cohort_unc[i], 1] +
                beta[cohort_unc[i], 2] * occupy_unc[i] +
                beta[cohort_unc[i], 3] * env_unc[i] +
                beta[cohort_unc[i], 4] * lit_unc[i] +
                beta[cohort_unc[i], 5] * (env_unc[i] * lit_unc[i]) +
                beta[cohort_unc[i], 6] * size_unc[i]) / alpha)));
    }
  }
  for(i in 1:N_cen) {
    increment_log_prob(weibull_ccdf_log(dur_cen[i], alpha,
          exp(-(beta[cohort_cen[i], 1] +
              beta[cohort_cen[i], 2] * occupy_cen[i] +
              beta[cohort_cen[i], 3] * env_cen[i] +
              beta[cohort_cen[i], 4] * lit_cen[i] +
              beta[cohort_unc[i], 5] * (env_cen[i] * lit_cen[i]) +
              beta[cohort_cen[i], 6] * size_cen[i]) / alpha)));
  }
}
generated quantities {
  vector[N] log_lik;

  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      log_lik[i] <- weibull_cdf_log(dur_unc[i], alpha, 
          exp(-(beta[cohort_unc[i], 1] +
              beta[cohort_unc[i], 2] * occupy_unc[i] +
              beta[cohort_unc[i], 3] * env_unc[i] +
              beta[cohort_unc[i], 4] * lit_unc[i] +
              beta[cohort_unc[i], 5] * (env_unc[i] * lit_unc[i]) +
              beta[cohort_unc[i], 6] * size_unc[i]) / alpha));
    } else {
      log_lik[i] <- weibull_log(dur_unc[i], alpha,
          exp(-(beta[cohort_unc[i], 1] +
              beta[cohort_unc[i], 2] * occupy_unc[i] +
              beta[cohort_unc[i], 3] * env_unc[i] +
              beta[cohort_unc[i], 4] * lit_unc[i] +
              beta[cohort_unc[i], 5] * (env_unc[i] * lit_unc[i]) +
              beta[cohort_unc[i], 6] * size_unc[i]) / alpha));
    }
  }
  for(i in 1:N_cen) {
    log_lik[i + N_unc] <- weibull_ccdf_log(dur_cen[i], alpha,
        exp(-(beta[cohort_cen[i], 1] +
            beta[cohort_cen[i], 2] * occupy_cen[i] +
            beta[cohort_cen[i], 3] * env_cen[i] +
            beta[cohort_cen[i], 4] * lit_cen[i] +
            beta[cohort_unc[i], 5] * (env_cen[i] * lit_cen[i]) +
            beta[cohort_cen[i], 6] * size_cen[i]) / alpha));
  }
}
