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

  int regime[O];
}
transformed data {
}
parameters {
  real<lower=0> alpha;
  real intercept;

  real group[C];
  real<lower=0> sigma_group;

  real cohort[O];
  real<lower=0> sigma_cohort[O];
  
  real cohort_mu[R];
  real<lower=0> sigma_regime;

  real<lower=0> sigma_cohcoh;
}
model {
  alpha ~ cauchy(0, 2.5);
  intercept ~ normal(0, 10);

  // class effect
  for(c in 1:C) {
    group[c] ~ normal(0, sigma_group);
  }
  sigma_group ~ cauchy(0, 2.5);

  // temporal effect
  for(o in 1:O) {
    cohort[o] ~ normal(cohort_mu[regime[o]], sigma_cohort[o]);
    sigma_cohort[o] ~ lognormal(0, sigma_cohcoh);
  }
  sigma_cohcoh ~ cauchy(0, 2.5);


  for(r in 1:R) {
    cohort_mu[r] ~ normal(0, sigma_regime);
  }
  sigma_regime ~ cauchy(0, 2.5);

  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      increment_log_prob(weibull_cdf_log(dur_unc[i], alpha, 
            exp(-(intercept + group[group_unc[i]] + cohort[cohort_unc[i]])/ alpha)));
    } else {
      increment_log_prob(weibull_log(dur_unc[i], alpha, 
            exp(-(intercept + group[group_unc[i]] + cohort[cohort_unc[i]])/ alpha)));
    }
  }

  for(i in 1:N_cen) {
    increment_log_prob(weibull_ccdf_log(dur_cen[i], alpha,
          exp(-(intercept + group[group_cen[i]] + cohort[cohort_cen[i]]) / alpha)));
  }
}
generated quantities {
}
