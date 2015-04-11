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
  real env_unc[N_unc];
  real<lower=0> dur_cen[N_cen];
  int group_cen[N_cen];
  int cohort_cen[N_cen];
  real env_cen[N_cen];

  real occupy_unc[N_unc];
  real occupy_cen[N_cen];
}
parameters {
  real<lower=0> alpha;
  vector[2] coef[O];
  vector[2] mu_prior;
  vector<lower=0>[2] sigma;
  corr_matrix[2] Omega;
}
transformed parameters {
  cov_matrix[2] Sigma;

  Sigma <- quad_form_diag(Omega, sigma);
}
model {
  alpha ~ cauchy(0, 2);

  // varying slopes, varying intercepts
  // done by cohort
  Omega ~ lkj_corr(2);
  sigma ~ cauchy(0, 1);
  for(i in 1:2) {
    mu_prior[i] ~ normal(0, 5);
  }
  for(i in 1:O) {
    coef[i] ~ multi_normal(mu_prior, Sigma);
  }

  // likelihood
  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      increment_log_prob(weibull_cdf_log(dur_unc[i], alpha, 
            exp(-(coef[cohort_unc[i], 1] +
                coef[cohort_unc[i], 2] * occupy_unc[i])/ alpha)));
    } else {
      increment_log_prob(weibull_log(dur_unc[i], alpha, 
            exp(-(coef[cohort_unc[i], 1] +
                coef[cohort_unc[i], 2] * occupy_unc[i])/ alpha)));
    }
  }
  for(i in 1:N_cen) {
    increment_log_prob(weibull_ccdf_log(dur_cen[i], alpha,
          exp(-(coef[cohort_cen[i], 1] +
              coef[cohort_cen[i], 2] * occupy_cen[i])/ alpha)));
  }
}
