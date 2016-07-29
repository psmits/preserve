functions {
  /**
   * return the log proability of observed fossil ate y given start, end, mid, and l
   *
   * @param y value between end and start of occurrence time
   * @param start time of origination
   * @param end time of loss
   * @param mid mode of the sampling distribution (default is (end - start) / 2)
   * @param l shape of the sampling distribution (default 4)
   * @return log probability of occurrence time
   */
  real pert_log(real y, real start, real end, real mid, real l) {
    real a1;
    real a2;
    real d;

    a1 = 1 + l * (mid - start) / (end - start);
    a2 = 1 + l * (end - mid) / (end - start);
    d = (y - start)^(a1 - 1) * (end - y)^(a2 - 1) / exp(lbeta(a1, a2)) / 
        (end - start)^(a1 + a2 - 1);

    return log(d);
  }
}
data {
  int N;  // total number of occurrences
  int S;  // total number of species
  int O;  // total number of cohorts

  real y[N];  // occurrence ages
  int sp[N];  // which species is occurrence n from
  vector[S] fad;  // species first appearance date
  vector[S] lad;  // species last appearance date

  int censored[N];
  
  real occupy[N];
  real env[N];
  real leng[N];

  vector[O] cohort_start;
  vector[O] cohort_end;
}
parameters {
  vector<lower=0,upper=1>[S] start_raw;
  vector<lower=1>[S] end_raw;

  real<lower=0> shape;

  // regression coefficients
  vector[5] mu_prior;
  vector[5] beta[O];  // cohort coef 
  corr_matrix[5] Omega;
  vector<lower=0>[5] sigma; 
}
transformed parameters {
  vector[S] start;  // actual origination time
  vector[S] end;  // actual exinction time
  vector[S] mid;  // mode
  
  cov_matrix[5] Sigma;
  
  // other misc transforms
  Sigma = quad_form_diag(Omega, sigma);

  // have to scale these because individual-level parameter constraints
  start = start_raw .* fad;
  end = end_raw .* lad;

  // constrain to be centered bell
  mid = ((end - start) / 2) + start;
}
model {
  vector[N] dur;
  vector[N] scale; // scale of weibull for end - start
  int cohort[S];
 
  // log absolute derivative of theta_raw .* d
  for(ss in 1:S) { // jacobian adjustment is just a constant
    target += (log(fabs(fad[ss])));
    target += (log(fabs(lad[ss])));
  }
  
  shape ~ lognormal(0, 0.3);

  // regression coefficients
  Omega ~ lkj_corr(2);
  sigma ~ cauchy(0, 1);
  mu_prior[1] ~ normal(0, 5);
  mu_prior[2] ~ normal(-1, 1);
  mu_prior[3] ~ normal(0, 1);
  mu_prior[4] ~ normal(1, 1);
  mu_prior[5] ~ normal(0, 1);

  for(i in 1:O) {
    beta[i] ~ multi_normal(mu_prior, Sigma);
  }
  
  for(i in 1:S) {
    scale[i] = exp(-(beta[cohort[i], 1] + 
          beta[cohort[i], 2] * occupy[i] + 
          beta[cohort[i], 3] * env[i] + 
          beta[cohort[i], 4] * (env[i]^2) +
          beta[cohort[i], 5] * leng[i])/ shape);
  }

  for(i in 1:S) {
    if((start[i] <= cohort_start[i]) && 
       (start[i] > cohort_end[i])) {
      cohort[i] = i;
    }
  }


  dur = end - start;
  // this is the meat of the whole operation
  target += weibull_lpdf(dur | shape, scale);


  for(n in 1:N) {
    y[n] ~ pert(start[sp[n]], end[sp[n]], mid[sp[n]], 4);
  }
}
