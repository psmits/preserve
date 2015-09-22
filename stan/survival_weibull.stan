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
  real alpha_trans[O];
  real alpha_mu;
  real<lower=0> alpha_sigma;

  vector[5] mu_prior;
  vector[5] beta[O];  // betas
  corr_matrix[5] Omega;
  vector<lower=0>[5] sigma;
  real<lower=0,upper=1> theta_species[N];
  real<lower=0,upper=1> theta_back[N];

  real<lower=0> tau[3];
  real<lower=0> lambda;
}
transformed parameters {
  real<lower=0> alpha[O];
  cov_matrix[5] Sigma;
  real env[N];

  for(i in 1:O) {
    alpha[i] <- exp(alpha_trans[i]);
  }

  Sigma <- quad_form_diag(Omega, sigma);

  for(i in 1:N) {
    env[i] <- theta_species[i] - theta_back[i];
  }
}
model {
  for(i in 1:O) {
    alpha_trans[i] ~ normal(alpha_mu, alpha_sigma);
  }
  alpha_mu ~ normal(0, 1);
  alpha_sigma ~ cauchy(0, 1);
  // varying slopes, varying intercepts
  // done by cohort
  Omega ~ lkj_corr(2);
  sigma ~ cauchy(0, 1);
  mu_prior[1] ~ normal(0, 5);
  mu_prior[2] ~ normal(-1, 1);
  mu_prior[3] ~ normal(0, tau[1] * lambda);
  mu_prior[4] ~ normal(0, tau[2] * lambda);
  mu_prior[5] ~ normal(0, tau[3] * lambda);
  tau ~ cauchy(0, 1);
  lambda ~ cauchy(0, 1);
  
  for(i in 1:O) {
    beta[i] ~ multi_normal(mu_prior, Sigma);
  }

  // uncertainty in environmental preference
  // here is the full posterior for each
  for(i in 1:N_unc) {
    theta_species[i] ~ beta(epi_unc[i] + 1, off_unc[i] + 1);
    theta_back[i] ~ beta(epi_bck_unc[i] + 1, off_bck_unc[i] + 1);
  }
  for(i in 1:N_cen) {
    theta_species[i + N_unc] ~ beta(epi_cen[i] + 1, off_cen[i] + 1);
    theta_back[i + N_unc] ~ beta(epi_bck_cen[i] + 1, off_bck_cen[i] + 1);
  }

  // likelihood
  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      increment_log_prob(weibull_cdf_log(dur_unc[i], alpha[cohort_unc[i]], 
            exp(-(beta[cohort_unc[i], 1] +
                beta[cohort_unc[i], 2] * occupy_unc[i] +
                beta[cohort_unc[i], 3] * env[i] +
                beta[cohort_unc[i], 4] * (env[i] * env[i]) +
                beta[cohort_unc[i], 5] * size_unc[i]) / alpha[cohort_unc[i]])));
    } else {
      increment_log_prob(weibull_log(dur_unc[i], alpha[cohort_unc[i]],
            exp(-(beta[cohort_unc[i], 1] +
                beta[cohort_unc[i], 2] * occupy_unc[i] +
                beta[cohort_unc[i], 3] * env[i] +
                beta[cohort_unc[i], 4] * (env[i] * env[i]) +
                beta[cohort_unc[i], 5] * size_unc[i]) / alpha[cohort_unc[i]])));
    }
  }
  for(i in 1:N_cen) {
    increment_log_prob(weibull_ccdf_log(dur_cen[i], alpha[cohort_cen[i]],
          exp(-(beta[cohort_cen[i], 1] +
              beta[cohort_cen[i], 2] * occupy_cen[i] +
              beta[cohort_cen[i], 3] * env[N_unc + i] +
              beta[cohort_cen[i], 4] * (env[N_unc + i] * env[N_unc + i]) +
              beta[cohort_cen[i], 5] * size_cen[i]) / alpha[cohort_cen[i]])));
  }
}
generated quantities {
  vector[N] log_lik;

  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      log_lik[i] <- weibull_cdf_log(dur_unc[i], alpha[cohort_unc[i]], 
          exp(-(beta[cohort_unc[i], 1] +
              beta[cohort_unc[i], 2] * occupy_unc[i] +
              beta[cohort_unc[i], 3] * env[i] +
              beta[cohort_unc[i], 4] * (env[i] * env[i]) +
              beta[cohort_unc[i], 5] * size_unc[i]) / alpha[cohort_unc[i]]));
    } else {
      log_lik[i] <- weibull_log(dur_unc[i], alpha[cohort_unc[i]],
          exp(-(beta[cohort_unc[i], 1] +
              beta[cohort_unc[i], 2] * occupy_unc[i] +
              beta[cohort_unc[i], 3] * env[i] +
              beta[cohort_unc[i], 4] * (env[i] * env[i]) +
              beta[cohort_unc[i], 5] * size_unc[i]) / alpha[cohort_unc[i]]));
    }
  }
  for(i in 1:N_cen) {
    log_lik[i + N_unc] <- weibull_ccdf_log(dur_cen[i], alpha[cohort_cen[i]],
        exp(-(beta[cohort_cen[i], 1] +
            beta[cohort_cen[i], 2] * occupy_cen[i] +
            beta[cohort_cen[i], 3] * env[N_unc + i] +
            beta[cohort_cen[i], 4] * (env[N_unc + i] * env[N_unc + i]) +
            beta[cohort_cen[i], 5] * size_cen[i]) / alpha[cohort_cen[i]]));
  }
}
