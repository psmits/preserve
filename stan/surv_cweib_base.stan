data {
  int N;  // number of samples
  int O;  // origination cohort

  real<lower=0> dur[N];  // duration

  int censored[N];  // 1 == censored, 0 == uncensored
 
  int cohort[N];
  real occupy[N];
  real env[N];
  real leng[N];

  real nocc[N];
}
parameters {
  real alpha_trans;

  // regression coefficients
  vector[5] mu_prior;
  vector[5] beta[O];  // cohort coef 
  corr_matrix[5] Omega;
  vector<lower=0>[5] sigma; 
  real delta;  // effect of abundance
}
transformed parameters {
  real<lower=0> alpha;
  cov_matrix[5] Sigma;

  vector[N] hold;
  
  // other misc transformations
  alpha = exp(alpha_trans);
  Sigma = quad_form_diag(Omega, sigma);

  // accumulate scale parameter
  for(i in 1:N) {
    hold[i] = exp(-(beta[cohort[i], 1] + 
          beta[cohort[i], 2] * occupy[i] + 
          beta[cohort[i], 3] * env[i] + 
          beta[cohort[i], 4] * (env[i]^2) +
          beta[cohort[i], 5] * leng[i] +
          delta * nocc[i])/ alpha);
  }
}
model {
  alpha_trans ~ normal(0, 0.5);

  // regression coefficients
  Omega ~ lkj_corr(2);
  sigma ~ cauchy(0, 1);
  mu_prior[1] ~ normal(0, 5);
  mu_prior[2] ~ normal(-1, 1);
  mu_prior[3] ~ normal(0, 1);
  mu_prior[4] ~ normal(1, 1);
  mu_prior[5] ~ normal(0, 0.5);

  delta ~ normal(0, 1);

  // this can be improved with non-centered parameterization
  //   adds parameters
  //   important of divergent transition issues
  for(i in 1:O) {
    beta[i] ~ multi_normal(mu_prior, Sigma);
  }

  // both estimate values of alp and bet (phi and lambda) for observed
  // then impute unobserved


  // likelihood / sampling statements
  for(i in 1:N) {
    if(censored[i] == 0) {
      if(dur[i] == 1) {
        target += weibull_lcdf(dur[i] | alpha, hold[i]);
      } else {
        target += weibull_lpdf(dur[i] | alpha, hold[i]);
      }
    } else {
      target += weibull_lccdf(dur[i] | alpha, hold[i]);
    }
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_tilde;

  // log_lik
  for(i in 1:N) {
    if(censored[i] == 0) {
      if(dur[i] == 1) {
        log_lik[i] = weibull_lcdf(dur[i] | alpha, hold[i]);
      } else {
        log_lik[i] = weibull_lpdf(dur[i] | alpha, hold[i]);
      }
    } else {
      log_lik[i] = weibull_lccdf(dur[i] | alpha, hold[i]);
    }
  }

  // posterior predictive simulations
  for(i in 1:N) {
    y_tilde[i] = weibull_rng(alpha, hold[i]);
  }
}
