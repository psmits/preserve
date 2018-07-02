data {
  int N;  // number of samples
  int O;  // origination cohort
  int N_obs;  // number observations will fully observed covariates
  int N_imp;  // number observations without fully observed covariates
  int K; // number of predictors
  int L;  // # group level predictors

  real<lower=0> dur[N];  // duration

  int censored[N];  // 1 == censored, 0 == uncensored
  int inclusion[N];  // 1 == fully observed, 0 == not observed
  int obs_ord[N_obs];
  int imp_ord[N_imp];
 
  int cohort[N];
  real occupy[N];
  real env[N];
  real leng[N];

  real<lower=0,upper=1> samp_obs[N_obs];  // observed sampling values
}
transformed data {
  matrix[O, L] u;

  for(o in 1:O) {
    for(l in 1:L) {
      u[o, l] = 1;
    }
  }
}
parameters {
  real alpha_trans;

  // regression coefficients
  // K = # covariates, O = # cohorts
  vector<lower=0>[K] tau;  // prior scale
  // no group-level predictors e.g. temp
  matrix[L, K] mu;
  matrix[K, O] z;
  cholesky_factor_corr[K] L_Omega;  // prior correlation

  // no prediction error scale because binary

  real delta;
  vector[4] gamma;
  real<lower=0,upper=1> samp_imp[N_imp];  // imputed sampling values
  real<lower=0.1> lambda;  // beta total count
}
transformed parameters {
  real<lower=0> alpha;

  real<lower=0,upper=1> phi[N];  // beta mean
  real<lower=0> alp[N];
  real<lower=0> bet[N];
  real<lower=0,upper=1> samp[N];
  
  vector[N] hold;
  
  matrix[O, K] beta;  // indiv K coef by group 0
  
  // other misc transformations
  alpha = exp(alpha_trans);

  // predict mean of inclusion function
  for(n in 1:N) {
    phi[n] = inv_logit(gamma[1] + 
        gamma[2] * occupy[n] + 
        gamma[3] * env[n] + 
        gamma[4] * leng[n]);
  }
  
  // parameter transformation
  for(n in 1:N) {
    alp[n] = lambda * phi[n];
    bet[n] = lambda * (1 - phi[n]);
  }
  
  // combining data and parameters
  samp[obs_ord] = samp_obs;
  samp[imp_ord] = samp_imp;
  
  // accumulate scale parameter
  for(i in 1:N) {
    hold[i] = exp(-(beta[cohort[i], 1] + 
          beta[cohort[i], 2] * occupy[i] + 
          beta[cohort[i], 3] * env[i] + 
          beta[cohort[i], 4] * (env[i]^2) +
          beta[cohort[i], 5] * leng[i] +
          delta * samp[i])/ alpha);
  }

  // coef madness
  beta = u * mu + (diag_pre_multiply(tau, L_Omega) * z)';
}
model {
  
  alpha_trans ~ normal(0, 1);

  // priors for group-level average effect
  tau ~ cauchy(0, 1);
  to_vector(z) ~ normal(0, 1);
  to_vector(mu) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(2);

  // both estimate values of alp and bet (phi and lambda) for observed
  // then impute unobserved
  lambda ~ pareto(0.1, 1.5);
  gamma ~ normal(0, 1);
  samp_obs ~ beta(alp[obs_ord], bet[obs_ord]);
  samp_imp ~ beta(alp[imp_ord], bet[imp_ord]);
  // this is creating a log(0) situation...


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
  matrix[K, K] Omega;
  matrix[K, K] Sigma;

  Omega = L_Omega * L_Omega';
  Sigma = quad_form_diag(Omega, tau);

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
    if(censored[i] == 1) {
      //y_tilde[i] = y_tilde[i] * uniform_rng(0, 1);
      if(y_tilde[i] > dur[i]) y_tilde[i] = dur[i];
    }
    if(dur[i] == 1 && censored[i] == 0 && y_tilde[i] < 1) {
      y_tilde[i] = 1;
    }
  }
}


