data {
  int N;  // number of samples
  int O;  // origination cohort
  int R;  // background regime
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

  int regime[O];
}
transformed data {
}
parameters {
  real<lower=0> alpha;
  real intercept;
  real slope;
  real<lower=0,upper=1> x_unc[N_unc];
  real<lower=0,upper=1> x_cen[N_cen];
}
model {
  alpha ~ cauchy(0, 2.5);
  intercept ~ normal(0, 10);
  slope ~ normal(0, 5);

  for(i in 1:N_unc) {
    x_unc ~ beta(1, 1);
    x_unc[i] ~ beta(epi_unc[i] + epi_bck_unc[i], off_unc[i] + off_bck_unc[i]);
  }
  for(i in 1:N_cen) {
    x_cen ~ beta(1, 1);
    x_cen[i] ~ beta(epi_cen[i] + epi_bck_cen[i], off_cen[i] + off_bck_cen[i]);
  }

  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      increment_log_prob(weibull_cdf_log(dur_unc[i], alpha, 
            exp(-(intercept + slope * x_unc[i])/ alpha)));
    } else {
      increment_log_prob(weibull_log(dur_unc[i], alpha, 
            exp(-(intercept + slope * x_unc[i])/ alpha)));
    }
  }

  for(i in 1:N_cen) {
    increment_log_prob(weibull_ccdf_log(dur_cen[i], alpha,
          exp(-(intercept + slope * x_cen[i])/ alpha)));
  }
}
