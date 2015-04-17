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
  real<lower=0> alpha;
  vector[4] mu_prior;
  vector[4] beta[O];  // betas
  corr_matrix[4] Omega;
  vector<lower=0>[4] sigma;
  vector<lower=0,upper=1>[N] pref;
}
transformed parameters {
  cov_matrix[4] Sigma;

  Sigma <- quad_form_diag(Omega, sigma);

  // uncertainty in environment
  // logit transform then rescale
}
model {
  real env[N];
  real pref_trans[N];
  for(i in 1:N) {
    pref_trans[i] <- logit(pref[i]);
  }
  for(i in 1:N) {
    env[i] <- (pref_trans[i] - mean(pref_trans)) / (2 * sd(pref_trans));
  }
  
  alpha ~ cauchy(0, 2);
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
      increment_log_prob(weibull_cdf_log(dur_unc[i], alpha, 
            exp(-(beta[cohort_unc[i], 1] +
                beta[cohort_unc[i], 2] * occupy_unc[i] +
                beta[cohort_unc[i], 3] * env_unc[i] +
                beta[cohort_unc[i], 4] * size_unc[i]) / alpha)));
    } else {
      increment_log_prob(weibull_log(dur_unc[i], alpha,
            exp(-(beta[cohort_unc[i], 1] +
                beta[cohort_unc[i], 2] * occupy_unc[i] +
                beta[cohort_unc[i], 3] * env_unc[i] +
                beta[cohort_unc[i], 4] * size_unc[i]) / alpha)));
    }
  }
  for(i in 1:N_cen) {
    increment_log_prob(weibull_ccdf_log(dur_cen[i], alpha,
          exp(-(beta[cohort_cen[i], 1] +
              beta[cohort_cen[i], 2] * occupy_cen[i] +
              beta[cohort_cen[i], 3] * env_cen[i] +
              beta[cohort_cen[i], 4] * size_cen[i]) / alpha)));
  }
}
