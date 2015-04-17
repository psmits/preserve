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
  int epi_unc[N_unc];
  int off_unc[N_unc];
  int epi_bck_unc[N_unc];
  int off_bck_unc[N_unc];
  real size_unc[N_unc];
  real<lower=0> dur_cen[N_cen];
  int group_cen[N_cen];
  int cohort_cen[N_cen];
  real occupy_cen[N_cen];
  real env_cen[N_cen];
  int epi_cen[N_cen];
  int off_cen[N_cen];
  int epi_bck_cen[N_cen];
  int off_bck_cen[N_cen];
  real size_cen[N_cen];
}
parameters {
  vector[4] mu_prior;
  vector[4] beta[O];  // betas
  corr_matrix[4] Omega;
  vector<lower=0>[4] sigma;
  vector<lower=0,upper=1>[N] pref;
}
transformed parameters {
  cov_matrix[4] Sigma;
  real pref_trans[N];

  Sigma <- quad_form_diag(Omega, sigma);
  
  for(i in 1:N) {
    pref_trans[i] <- logit(pref[i]);
  }

  // uncertainty in environment
  // logit transform then rescale
}
model {
  real env[N];
  for(i in 1:N) {
    env[i] <- (pref_trans[i] - mean(pref_trans)) / (2 * sd(pref_trans));
  }
  // varying slopes, varying intercepts
  // done by cohort
  Omega ~ lkj_corr(2);
  sigma ~ cauchy(0, 1);
  for(i in 1:4) {
    mu_prior[i] ~ normal(0, 5);
  }
  for(i in 1:O) {
    beta[i] ~ multi_normal(mu_prior, Sigma);
  }

# uncertainty in environmental preference
# here is the full posterior
  for(n in 1:N_unc) {
    pref[n] ~ beta(epi_unc[n] + epi_bck_unc[n], off_unc[n] + off_bck_unc[n]);
  }
  for(n in 1:N_cen) {
    pref[n + N_unc] ~ beta(epi_cen[n] + epi_bck_cen[n], off_cen[n] + off_bck_cen[n]);
  }

  // likelihood
  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      increment_log_prob(exponential_cdf_log(dur_unc[i],
            exp(beta[cohort_unc[i], 1] +
              beta[cohort_unc[i], 2] * occupy_unc[i] +
              beta[cohort_unc[i], 3] * pref_trans[i] +
              beta[cohort_unc[i], 4] * size_unc[i])));
    } else {
      increment_log_prob(exponential_log(dur_unc[i],
            exp(beta[cohort_unc[i], 1] +
              beta[cohort_unc[i], 2] * occupy_unc[i] +
              beta[cohort_unc[i], 3] * pref_trans[i] +
              beta[cohort_unc[i], 4] * size_unc[i])));
    }
  }
  for(i in 1:N_cen) {
    increment_log_prob(exponential_ccdf_log(dur_cen[i],
          exp(beta[cohort_cen[i], 1] +
            beta[cohort_cen[i], 2] * occupy_cen[i] +
            beta[cohort_cen[i], 3] * pref_trans[N_unc + i] +
            beta[cohort_cen[i], 4] * size_cen[i])));
  }
}
generated quantities {
  vector[N] log_lik;

  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      log_lik[i] <- exponential_cdf_log(dur_unc[i],
          exp(beta[cohort_unc[i], 1] +
            beta[cohort_unc[i], 2] * occupy_unc[i] +
            beta[cohort_unc[i], 3] * pref_trans[i] +
            beta[cohort_unc[i], 4] * size_unc[i]));
    } else {
      log_lik[i] <- exponential_log(dur_unc[i],
          exp(beta[cohort_unc[i], 1] +
            beta[cohort_unc[i], 2] * occupy_unc[i] +
            beta[cohort_unc[i], 3] * pref_trans[i] +
            beta[cohort_unc[i], 4] * size_unc[i]));
    }
  }
  for(i in 1:N_cen) {
    log_lik[i + N_unc] <- exponential_ccdf_log(dur_cen[i],
        exp(beta[cohort_cen[i], 1] +
          beta[cohort_cen[i], 2] * occupy_cen[i] +
          beta[cohort_cen[i], 3] * pref_trans[N_unc + i] +
          beta[cohort_cen[i], 4] * size_cen[i]));
  }
}
