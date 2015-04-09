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
  real<lower=0> dur_cen[N_cen];
  int group_cen[N_cen];
  int cohort_cen[N_cen];

  int epi_unc[N_unc];
  int off_unc[N_unc];
  int epi_bck_unc[N_unc];
  int off_bck_unc[N_unc];
  int epi_cen[N_cen];
  int off_cen[N_cen];
  int epi_bck_cen[N_cen];
  int off_bck_cen[N_cen];
}
transformed data {
  vector[2] mu_prior;
  for(i in 1:2) mu_prior[i] <- 0;
}
parameters {
  real<lower=0> alpha;

  vector[2] coef[C];
  vector<lower=0>[2] sigma;
  corr_matrix[2] Omega;

  real temp[O];
  real<lower=0> tau;

  real<lower=0,upper=1> x_unc[N_unc];
  real<lower=0,upper=1> x_cen[N_cen];
}
transformed parameters {
  cov_matrix[2] Sigma;
  Sigma <- quad_form_diag(Omega, sigma);
}
model {
  alpha ~ cauchy(0, 2.5);
  
  // varying slopes, varying intercepts
  Omega ~ lkj_corr(2);
  sigma ~ cauchy(0, 1);
  for(i in 1:C) {
    coef[i] ~ multi_normal(mu_prior, Sigma);
  }

  // non-nested temporal effect
  for(i in 1:O) {
    temp[i] ~ normal(0, tau);
  }
  tau ~ cauchy(0, 1);
  
  // posterior of environment preference
  for(i in 1:N_unc) {
    x_unc[i] ~ beta(epi_unc[i] + epi_bck_unc[i], off_unc[i] + off_bck_unc[i]);
  }
  for(i in 1:N_cen) {
    x_cen[i] ~ beta(epi_cen[i] + epi_bck_cen[i], off_cen[i] + off_bck_cen[i]);
  }

  // likelihood
  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      increment_log_prob(weibull_cdf_log(dur_unc[i], alpha, 
            exp(-(coef[group_unc[i], 1] + temp[cohort_unc[i]] +
                coef[group_unc[i], 2] * x_unc[i])/ alpha)));
    } else {
      increment_log_prob(weibull_log(dur_unc[i], alpha, 
            exp(-(coef[group_unc[i], 1] + temp[cohort_unc[i]] +
                coef[group_unc[i], 2] * x_unc[i])/ alpha)));
    }
  }
  for(i in 1:N_cen) {
    increment_log_prob(weibull_ccdf_log(dur_cen[i], alpha,
          exp(-(coef[group_cen[i], 1] + temp[cohort_cen[i]] +
              coef[group_cen[i], 2] * x_cen[i])/ alpha)));
  }
}
