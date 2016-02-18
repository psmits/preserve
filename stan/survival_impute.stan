data {
  int N;  // number of samples
  int O;  // origination cohort
  int N_obs;
  int N_imp;
  
  real<lower=0> dur[N];  // duration of uncensored obs

  int censored[N];
  int inclusion[N];
 
  int cohort[N];
  real occupy[N];
  real env[N];
  real leng[N];

  real samp_obs[N_obs];
}
parameters {
  real alpha_trans;

  // regression coefficients
  vector[6] mu_prior;
  vector[6] beta[O];  // cohort coef 
  corr_matrix[6] Omega;
  vector<lower=0>[6] sigma; 

  vector[4] gamma;
  real samp_imp[N_imp];
  real<lower=0.1> lambda;
}
transformed parameters {
  real<lower=0> alpha;  // individual shape
  cov_matrix[6] Sigma;

  real<lower=0,upper=1> phi[N];
  real alp[N];
  real bet[N];
  
  real samp[N];
  
  // parameter transformation
  for(n in 1:N) {
    alp[n] <- lambda * phi[n];
    bet[n] <- lambda * (1 - phi[n]);
  }
  
  // other misc transformations
  alpha <- exp(alpha_trans);
  Sigma <- quad_form_diag(Omega, sigma);

  for(n in 1:N) {
    if(inclusion[n] == 1) {
      samp[n] <- samp_obs[n];
    } else {
      samp[n] <- samp_imp[n];
    }
  }

  for(n in 1:N) {
    phi[n] <- gamma[1] + 
      gamma[2] * occupy[n] + 
      gamma[3] * env[n] + 
      gamma[4] * leng[n];
  }
}
model {
  
  alpha_trans ~ normal(0, 1);

  // regression coefficients
  Omega ~ lkj_corr(2);
  sigma ~ cauchy(0, 1);
  mu_prior[1] ~ normal(0, 5);
  mu_prior[2] ~ normal(-1, 1);
  mu_prior[3] ~ normal(0, 1);
  mu_prior[4] ~ normal(1, 1);
  mu_prior[5] ~ normal(0, 1);
  mu_prior[6] ~ normal(0, 1);

  for(i in 1:O) {
    beta[i] ~ multi_normal(mu_prior, Sigma);
  }

  // have to update phi to be a regression
  lambda ~ pareto(0.1, 1.5);
  gamma ~ normal(0, 1);
  for(n in 1:N) {
    if(inclusion[n] == 1) {
      samp_obs[n] ~ beta(alp[n], bet[n]);
    } else {
      samp_imp[n] ~ beta(alp[n], bet[n]);
    }
  }

  // likelihood / sampling statements
  for(i in 1:N) {
    if(censored[i] == 0) {
      if(dur[i] == 1) {
        increment_log_prob(weibull_cdf_log(dur[i], alpha,
              exp(-(beta[cohort[i], 1] + 
                  beta[cohort[i], 2] * occupy[i] + 
                  beta[cohort[i], 3] * env[i] + 
                  beta[cohort[i], 4] * (env[i]^2) +
                  beta[cohort[i], 5] * leng[i] + 
                  beta[cohort[i], 6] * samp[i])/ alpha)));
      } else {
        increment_log_prob(weibull_log(dur[i], alpha,
              exp(-(beta[cohort[i], 1] + 
                  beta[cohort[i], 2] * occupy[i] + 
                  beta[cohort[i], 3] * env[i] + 
                  beta[cohort[i], 4] * (env[i]^2) +
                  beta[cohort[i], 5] * leng[i] + 
                  beta[cohort[i], 6] * samp[i])/ alpha)));
      }
    } else {
      increment_log_prob(weibull_ccdf_log(dur[i], alpha,
            exp(-(beta[cohort[i], 1] + 
                beta[cohort[i], 2] * occupy[i] + 
                beta[cohort[i], 3] * env[i] + 
                beta[cohort[i], 4] * (env[i]^2) +
                beta[cohort[i], 5] * leng[i] + 
                beta[cohort[i], 6] * samp[i])/ alpha)));
    }
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_tilde;
  vector[N] hold;

  for(i in 1:N) {
    hold[i] <- exp(-(beta[cohort[i], 1] +
          beta[cohort[i], 2] * occupy[i] +
          beta[cohort[i], 3] * env[i] + 
          beta[cohort[i], 4] * (env[i]^2) +
          beta[cohort[i], 5] * leng[i] +
          beta[cohort[i], 6] * samp[i]) / alpha);
  }

  // log_lik
  for(i in 1:N) {
    if(censored[i] == 0) {
      if(dur[i] == 1) {
        log_lik[i] <- weibull_cdf_log(dur[i], alpha, hold[i]);
      } else {
        log_lik[i] <- weibull_log(dur[i], alpha, hold[i]);
      }
    } else {
      log_lik[i] <- weibull_ccdf_log(dur[i], 
          alpha, hold[i]);
    }
  }

  // posterior predictive simulations
  for(i in 1:N) {
    y_tilde[i] <- weibull_rng(alpha, hold[i]);
  }
}
