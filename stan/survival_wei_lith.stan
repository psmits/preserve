data {
  int N;  // # taxa
  int O;  // # cohort
  int N_unc;  // # uncensored
  int N_cen;  // # censored

  int cohort_unc[N_unc];  // cohort membership
  int cohort_cen[N_cen];

  real<lower=0> dur_unc[N_unc];  // duration
  real<lower=0> dur_cen[N_cen];

  real occupy_unc[N_unc];  // range size
  real occupy_cen[N_cen];

  real size_unc[N_unc];  // body size
  real size_cen[N_cen];

  int epi_unc[N_unc];  // environment; epi = epicontinental, off = open ocean
  int epi_cen[N_cen];
  int off_unc[N_unc];
  int off_cen[N_cen];
  int epi_bck_unc[N_unc];
  int epi_bck_cen[N_cen];
  int off_bck_unc[N_unc];
  int off_bck_cen[N_cen];

  int car_unc[N_unc];  // lithology; car = carbonate, cla = clastic
  int car_cen[N_cen];
  int cla_unc[N_unc];
  int cla_cen[N_cen];
  int car_bck_unc[N_unc];
  int car_bck_cen[N_cen];
  int cla_bck_unc[N_unc];
  int cla_bck_cen[N_cen];
}
parameters {
  real<lower=0> alpha;
  vector[9] mu_prior;
  vector[9] beta[O];  // betas
  corr_matrix[9] Omega;
  vector<lower=0>[9] sigma;
  real<lower=0,upper=1> theta_species[N];
  real<lower=0,upper=1> theta_back[N];
  real<lower=0,upper=1> phi_species[N];
  real<lower=0,upper=1> phi_back[N];

  real<lower=0> lambda;
  real<lower=0> tau[6];
}
transformed parameters {
  cov_matrix[9] Sigma;
  real env[N];
  real lit[N];

  Sigma <- quad_form_diag(Omega, sigma);

  for(i in 1:N) {
    env[i] <- theta_species[i] - theta_back[i];
    lit[i] <- phi_species[i] - phi_back[i];
  }
}
model {
  alpha ~ cauchy(0, 2);
  // varying slopes, varying intercepts
  // done by cohort
  Omega ~ lkj_corr(2);
  sigma ~ cauchy(0, 1);
  for(i in 1:9) {
    if(i == 2) {
      mu_prior[i] ~ normal(-1, 1);
    } else if(i > 3) {
      mu_prior[i] ~ normal(0, tau[i - 3] * lambda);
    } else {
      mu_prior[i] ~ normal(0, 5);
    }
  }
  tau ~ cauchy(0, 1);
  lambda ~ cauchy(0, 1);
  
  for(i in 1:O) {
    beta[i] ~ multi_normal(mu_prior, Sigma);
  }

  // uncertainty in environmental preference
  // here is the full posterior
  for(i in 1:N_unc) {
    theta_species[i] ~ beta(epi_unc[i] + 1, off_unc[i] + 1);
    theta_back[i] ~ beta(epi_bck_unc[i] + 1, off_bck_unc[i] + 1);
    phi_species[i] ~ beta(car_unc[i] + 1, cla_unc[i] + 1);
    phi_back[i] ~ beta(car_bck_unc[i] + 1, cla_bck_unc[i] + 1);
  }
  for(i in 1:N_cen) {
    theta_species[i + N_unc] ~ beta(epi_cen[i] + 1, off_cen[i] + 1);
    theta_back[i + N_unc] ~ beta(epi_bck_cen[i] + 1, off_bck_cen[i] + 1);
    phi_species[i + N_unc] ~ beta(car_cen[i] + 1, cla_cen[i] + 1);
    phi_back[i + N_unc] ~ beta(car_bck_cen[i] + 1, cla_bck_cen[i] + 1);
  }

  // likelihood
  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      increment_log_prob(weibull_cdf_log(dur_unc[i], alpha,
            exp(-(beta[cohort_unc[i], 1] +
              beta[cohort_unc[i], 2] * occupy_unc[i] +
              beta[cohort_unc[i], 3] * size_unc[i] +
              beta[cohort_unc[i], 4] * env[i] +
              beta[cohort_unc[i], 5] * env[i]^2 +
              beta[cohort_unc[i], 6] * lit[i] +
              beta[cohort_unc[i], 7] * lit[i]^2 +
              beta[cohort_unc[i], 8] * (env[i] * lit[i]) +
              beta[cohort_unc[i], 9] * (env[i] * lit[i])^2) / alpha)));
    } else {
      increment_log_prob(weibull_log(dur_unc[i], alpha,
            exp(-(beta[cohort_unc[i], 1] +
              beta[cohort_unc[i], 2] * occupy_unc[i] +
              beta[cohort_unc[i], 3] * size_unc[i] +
              beta[cohort_unc[i], 4] * env[i] +
              beta[cohort_unc[i], 5] * env[i]^2 +
              beta[cohort_unc[i], 6] * lit[i] +
              beta[cohort_unc[i], 7] * lit[i]^2 +
              beta[cohort_unc[i], 8] * (env[i] * lit[i]) +
              beta[cohort_unc[i], 9] * (env[i] * lit[i])^2) 
              / alpha)));
    }
  }
  for(i in 1:N_cen) {
    increment_log_prob(weibull_ccdf_log(dur_cen[i], alpha,
          exp(-(beta[cohort_cen[i], 1] +
            beta[cohort_cen[i], 2] * occupy_cen[i] +
            beta[cohort_cen[i], 3] * size_cen[i] +
            beta[cohort_cen[i], 4] * env[N_unc + i] +
            beta[cohort_cen[i], 5] * env[N_unc + i]^2 +
            beta[cohort_cen[i], 6] * lit[N_unc + i] +
            beta[cohort_cen[i], 7] * lit[N_unc + i]^2 +
            beta[cohort_cen[i], 8] * (env[N_unc + i] * lit[N_unc + i]) +
            beta[cohort_cen[i], 9] * (env[N_unc + i] * lit[N_unc + i])^2) 
            / alpha)));
  }
}
generated quantities {
  vector[N] log_lik;

  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      log_lik[i] <- weibull_cdf_log(dur_unc[i], alpha,
          exp(-(beta[cohort_unc[i], 1] +
            beta[cohort_unc[i], 2] * occupy_unc[i] +
            beta[cohort_unc[i], 3] * size_unc[i] +
            beta[cohort_unc[i], 4] * env[i] +
            beta[cohort_unc[i], 5] * env[i]^2 +
            beta[cohort_unc[i], 6] * lit[i] +
            beta[cohort_unc[i], 7] * lit[i]^2 +
            beta[cohort_unc[i], 8] * (env[i] * lit[i]) +
            beta[cohort_unc[i], 9] * (env[i] * lit[i])^2) 
            / alpha));
    } else {
      log_lik[i] <- weibull_log(dur_unc[i], alpha,
          exp(-(beta[cohort_unc[i], 1] +
            beta[cohort_unc[i], 2] * occupy_unc[i] +
            beta[cohort_unc[i], 3] * size_unc[i] +
            beta[cohort_unc[i], 4] * env[i] +
            beta[cohort_unc[i], 5] * env[i]^2 +
            beta[cohort_unc[i], 6] * lit[i] +
            beta[cohort_unc[i], 7] * lit[i]^2 +
            beta[cohort_unc[i], 8] * (env[i] * lit[i]) +
            beta[cohort_unc[i], 9] * (env[i] * lit[i])^2) 
            / alpha));
    }
  }
  for(i in 1:N_cen) {
    log_lik[i + N_unc] <- weibull_ccdf_log(dur_cen[i], alpha,
        exp(-(beta[cohort_cen[i], 1] +
          beta[cohort_cen[i], 2] * occupy_cen[i] +
          beta[cohort_cen[i], 3] * size_cen[i] +
          beta[cohort_cen[i], 4] * env[N_unc + i] +
          beta[cohort_cen[i], 5] * env[N_unc + i]^2 +
          beta[cohort_cen[i], 6] * lit[N_unc + i] +
          beta[cohort_cen[i], 7] * lit[N_unc + i]^2 +
          beta[cohort_cen[i], 8] * (env[N_unc + i] * lit[N_unc + i]) +
          beta[cohort_cen[i], 9] * (env[N_unc + i] * lit[N_unc + i])^2) 
          / alpha));
  }
}

