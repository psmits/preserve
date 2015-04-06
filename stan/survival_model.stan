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
  real cohort_mu[R];
  real<lower=0> sigma_cohort;
  real<lower=0> sigma_regime;
}
model {
  alpha ~ cauchy(0, 2.5);
  intercept ~ normal(0, 10);

  // class effect
  for(c in 1:C) {
    group[c] ~ normal(0, sigma_group);
  }
  sigma_group ~ cauchy(0, 1);

  // temporal effect
  for(o in 1:O) {
    cohort[o] ~ normal(cohort_mu[regime[o]], sigma_cohort);
  }
  sigma_cohort ~ cauchy(0, 1);

  for(r in 1:R) {
    cohort_mu[r] ~ normal(0, sigma_regime);
  }
  sigma_regime ~ cauchy(0, 1);

  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      increment_log_prob(weibull_cdf_log(dur_unc[i], alpha, 
            exp(-(intercept + 
            group[group_unc[i]] + cohort[cohort_unc[i]])/ alpha)));
    } else {
      increment_log_prob(weibull_log(dur_unc[i], alpha, 
            exp(-(intercept + 
            group[group_unc[i]] + cohort[cohort_unc[i]])/ alpha)));
    }
  }

  for(i in 1:N_cen) {
    increment_log_prob(weibull_ccdf_log(dur_cen[i], alpha,
          exp(-(intercept + 
            group[group_cen[i]] + cohort[cohort_cen[i]])/ alpha)));
  }
}
generated quantities {
  vector[N] log_lik;

  for(i in 1:N_unc) {
    if(dur_unc[i] == 1) {
      log_lik[i] <- weibull_cdf_log(dur_unc[i], alpha,
          exp(-(intercept + 
              group[group_unc[i]] + cohort[cohort_unc[i]])/ alpha));
    } else {
      log_lik[i] <- weibull_log(dur_unc[i], alpha,
          exp(-(intercept + 
              group[group_unc[i]] + cohort[cohort_unc[i]])/ alpha));
    }
  }
  for(j in 1:N_cen) {
    log_lik[N_unc + j] <- weibull_ccdf_log(dur_cen[j], alpha,
        exp(-(intercept + 
              group[group_cen[j]] + cohort[cohort_cen[j]])/ alpha));
  }
}
