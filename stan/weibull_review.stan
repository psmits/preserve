data {
  int N;  // number of samples
  int O;  // origination cohort
  int N_unc;  // number uncensored obs
  int N_cen;  // number censored obs
  real<lower=0> dur_unc[N_unc];  // duration of uncensored obs
  int cohort_unc[N_unc];  // cohort membership
  real occupy_unc[N_unc];  // range size
  real size_unc[N_unc];  // body size
  int epi_unc[N_unc];  // number of epicontinental occ
  int off_unc[N_unc];  // number of open ocean occ
  real<lower=0> dur_cen[N_cen];
  int cohort_cen[N_cen];
  real occupy_cen[N_cen];
  int epi_cen[N_cen];
  int off_cen[N_cen];
  real size_cen[N_cen];
}
parameters {
  // measures environmental preference
  real<lower=0,upper=1> theta_taxon[N_unc + N_cen];
  real theta_mu;
  real<lower=0> sigma_theta;

  real<lower=0> alpha_trans[O];  // cohort shape
  real alpha_mu;  // shared (transformed) shape
  real<lower=0> sigma_alpha;

  // regression coefficients
  vector[5] mu_prior;
  vector[5] beta[O];  // cohort coef 
  corr_matrix[5] Omega;
  vector<lower=0>[5] sigma; 
}
transformed parameters {
  real theta_taxon_real[N_unc + N_cen];
  real<lower=0> alpha[O];
  cov_matrix[5] Sigma;

  for(n in 1:N) {
    theta_taxon_real[n] <- logit(theta_taxon[n]);
  }
  for(o in 1:O) {
    alpha[o] <- exp(alpha_trans[o]);
  }
  Sigma <- quad_form_diag(Omega, sigma);
}
model {
  // environmental preference
  theta_mu ~ normal(0, 1);
  sigma_theta ~ cauchy(0, 1);

  increment_log_prob(normal_log(theta_taxon_real, theta_mu, sigma_theta));
  for(n in 1:N_unc) {
    increment_log_prob(binomial_log(epi_unc[n], 
          epi_unc[n] + off_unc[n], theta_taxon[n]));
  }
  for(n in 1:N_cen) {
    increment_log_prob(binomial_log(epi_cen[n], 
          epi_cen[n] + off_cen[n], theta_taxon[N_unc + n]));
  }

  // hierarchical shape parameter
  alpha_mu ~ normal(0, 1);
  sigma_alpha ~ cauchy(0, 1);

  for(i in 1:O) {
    alpha_trans[i] ~ normal(alpha_mu, sigma_alpha);
  }

  // regression coefficients
  Omega ~ lkj_corr(2);
  sigma ~ cauchy(0, 1);
  mu_prior[1] ~ normal(0, 5);
  mu_prior[2] ~ normal(-1, 1);
  mu_prior[3] ~ normal(0, 1);
  mu_prior[4] ~ normal(0, 1);
  mu_prior[5] ~ normal(0, 1);

  for(i in 1:O) {
    beta[i] ~ multi_normal(mu_prior, Sigma);
  }

  // likelihood / sampling statements
  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      increment_log_prob(weibull_cdf_log(dur_unc[i], alpha[cohort_unc[i]],
            exp(-(beta[cohort_unc[i], 1] + 
                beta[cohort_unc[i], 2] * occupy_unc[i] + 
                beta[cohort_unc[i], 3] * theta_taxon[i] + 
                beta[cohort_unc[i], 4] + (theta_taxon[i] * theta_taxon[i]) + 
                beta[cohort_unc[i], 5] * size_unc[i]) / alpha[cohort_unc[i]])));
    } else {
      increment_log_prob(weibull_log(dur_unc[i], alpha[cohort_unc[i]],
            exp(-(beta[cohort_unc[i], 1] + 
                beta[cohort_unc[i], 2] * occupy_unc[i] + 
                beta[cohort_unc[i], 3] * theta_taxon[i] + 
                beta[cohort_unc[i], 4] + (theta_taxon[i] * theta_taxon[i]) + 
                beta[cohort_unc[i], 5] * size_unc[i]) / alpha[cohort_unc[i]])));
    }
  }
  for(i in 1:N_cen) {
    increment_log_prob(weibull_ccdf_log(dur_unc[i], alpha[cohort_unc[i]],
          exp(-(beta[cohort_cen[i], 1] + 
              beta[cohort_cen[i], 2] * occupy_cen[i] + 
              beta[cohort_cen[i], 3] * theta_taxon[i] + 
              beta[cohort_cen[i], 4] + (theta_taxon[i] * theta_taxon[i]) + 
              beta[cohort_cen[i], 5] * size_cen[i]) / alpha[cohort_cen[i]])));
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_tilde;
  vector[N] hold;

  for(i in 1:N_unc) {
    hold[i] <- exp(-(beta[cohort_unc[i], 1] +
          beta[cohort_unc[i], 2] * occupy_unc[i] +
          beta[cohort_unc[i], 3] * theta_taxon[i] + 
          beta[cohort_unc[i], 4] + (theta_taxon[i] * theta_taxon[i]) + 
          beta[cohort_unc[i], 5] * size_unc[i]) / alpha[cohort_unc[i]]);
  }
  for(i in 1:N_cen) {
    hold[i + N_unc] <- exp(-(beta[cohort_cen[i], 1] +
          beta[cohort_cen[i], 2] * occupy_cen[i] +
          beta[cohort_cen[i], 3] * theta_taxon[i] + 
          beta[cohort_cen[i], 4] + (theta_taxon[i] * theta_taxon[i]) + 
          beta[cohort_cen[i], 5] * size_unc[i]) / alpha[cohort_unc[i]]);
  }

  // log_lik
  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      log_lik[i] <- weibull_cdf_log(dur_unc[i], alpha[cohort_unc[i]], hold[i]);
    } else {
      log_lik[i] <- weibull_log(dur_unc[i], alpha[cohort_unc[i]], hold[i]);
    }
  }
  for(i in 1:N_cen) {
    log_lik[i + N_unc] <- weibull_ccdf_log(dur_cen[i], 
        alpha[cohort_cen[i]], hold[i]);
  }

  // posterior predictive simulations
  for(i in 1:N_unc) {
    y_tilde[i] <- weibull_rng(alpha[cohort_unc[i]], hold[i]);
  }
  for(i in 1:N_cen) {
    y_tilde[i + N_unc] <- weibull_rng(alpha[cohort_cen[i]], hold[i]);
  }
}
