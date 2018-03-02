functions {
  /**
   * Discrete Weibull distribution log probability mass function
   * 
   * @param y discrete value >= 0
   * @param alpha scale parameter > 0
   * @param beta shape parameter > 0
   * @return log probability of discrete time
   */
  real discrete_weibull_lpmf(int y, real alpha, real beta) {
    real o;

    o = exp(-(y / alpha) ^ beta) - exp(-((y + 1) / alpha) ^ beta);

    return log(o);
  }

  /**
   * Discrete Weibull distribution log cummulative distribution function
   * 
   * @param y discrete value >=0
   * @param alpha scale parameter > 0
   * @param beta shape parameter > 0
   * @return log CDF
   */
  real discrete_weibull_lcdf_lp(int y, real alpha, real beta) {
    real o;

    o = 1 - exp(-((y + 1) / alpha) ^ beta);

    return log(o);
  }

  /**
   * Discrete Weibull distribution log complementary cummulative distribution
   * function
   * 
   * @param y discrete value >= 0
   * @param alpha scale parameter > 0
   * @param beta shape parameter > 0
   * @return log CCDF
   */
  real discrete_weibull_lccdf_lp(int y, real alpha, real beta) {
    real o;

    o = 1 - exp(-((y + 1) / alpha) ^ beta);

    return log(1 - o);
  }

  /**
   * Survival function for discrete Weibull distribution
   * 
   * @param y discrete value >= 0
   * @param alpha scale parameter > 0
   * @param beta scale parameter > 0
   * @return probability of observed existing past y
   */
  real discrete_weibull_survival(int y, real alpha, real beta) {
    real o;
    real q;

    q = exp(-(alpha) ^ (-(beta)));

    o = q ^ (y  ^ beta);

    return o;
  }

  /**
   * Faliure "rate" function for discrete Weibull distribution
   * 
   * @param y discrete value >= 0
   * @param alpha scale parameter > 0
   * @param beta scale parameter > 0
   * @return probability of observed existing past y
   */
  real discrete_weibull_hazard(int y, real alpha, real beta) {
    real o;
    real q;

    q = exp(-(alpha) ^ (-(beta)));

    o = 1 - (q ^ (((y + 1) ^ beta) - (y ^ beta)));

    return o;
  }

  /**
   * Quantile function for discrete Weibull distribution
   * 
   * @param p real value 0 <= p <= 1
   * @param alpha scale parameter > 0
   * @param beta scale parameter > 0
   * @return quantile value from discrete weibull distribution 
   */
  real discrete_weibull_quantile(real p, real alpha, real beta) {
    real o;
    real q;

    q = exp(-(alpha) ^ (-(beta)));

    o = ceil((log(1 - p) / log(q)) ^ (1 / beta));

    return o;
  }

  /**
   * RNG function for discrete Weibull distribution
   * 
   * @param alpha scale parameter > 0
   * @param beta scale parameter > 0
   * @return random value from discrete weibull distribution 
   */
  real discrete_weibull_rng(real alpha, real beta) {
    real o;
    real u;

    u = uniform_rng(0, 1);

    o = discrete_weibull_quantile(u, alpha, beta);

    return o;
  }
}
data {
  int N;  // number of samples
  int O;  // origination cohort
  int N_obs;  // number observations will fully observed covariates
  int N_imp;  // number observations without fully observed covariates

  int<lower=0> dur[N];  // response: duration

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
parameters {

  // regression coefficients
  vector[5] mu_prior;
  vector[5] beta[O];  // cohort coef 
  corr_matrix[5] Omega;
  vector<lower=0>[5] sigma; 

  real delta;
  vector[4] gamma;
  real<lower=0,upper=1> samp_imp[N_imp];  // imputed sampling values
  real<lower=0.1> lambda;  // beta total count
}
transformed parameters {
  cov_matrix[5] Sigma;

  real<lower=0,upper=1> phi[N];  // beta mean
  real<lower=0> alp[N];
  real<lower=0> bet[N];
  real<lower=0,upper=1> samp[N];
  
  vector[N] hold;  // actual scale param for weibull
  
  // other misc transformations
  Sigma = quad_form_diag(Omega, sigma);

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
          delta * samp[i]));
  }
}
model {
  
  // regression coefficients
  Omega ~ lkj_corr(2);
  sigma ~ cauchy(0, 1);
  mu_prior[1] ~ normal(0, 5);
  mu_prior[2] ~ normal(-1, 1);
  mu_prior[3] ~ normal(0, 1);
  mu_prior[4] ~ normal(1, 1);
  mu_prior[5] ~ normal(0, 0.5);
  delta ~ normal(0, 1);

  for(i in 1:O) {
    beta[i] ~ multi_normal(mu_prior, Sigma);
  }

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
      target += discrete_weibull_lpmf(dur[i] | 1, hold[i]);
    } else {
      target += discrete_weibull_lccdf_lp(dur[i], 1, hold[i]);
    }
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_tilde;
  vector[N] survival_est;
  vector[N] hazard_est;

  for(i in 1:N) {
    // survival probability and hazard "rate"
    survival_est[i] = discrete_weibull_survival(dur[i], 1, hold[i]);
    hazard_est[i] = discrete_weibull_hazard(dur[i], 1, hold[i]);

    // log-likelihood calculation
    {
      real o;
      if(censored[i] == 0) {
        log_lik[i] = discrete_weibull_lpmf(dur[i] | 1, hold[i]);
      } else {
        o = 1 - exp(-((dur[i] + 1) / 1) ^ hold[i]);
        log_lik[i] = log(1 - o);
      }
    }

    // posterior predictive simulations
    y_tilde[i] = discrete_weibull_rng(1, hold[i]);
  }
}

