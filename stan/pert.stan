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
    real a;
    real b;
    real d;

    a = ((mid - end) * l) / (start - end);
    b = ((start - mid) * l) / (start - end);
    d = (((start - end)^(-1 - l)) * ((start - y)^a) * (-end + y)^b) /
      ((tgamma(a + 1) * tgamma(b + 1)) / (tgamma(a + b + 2)));

    return log(d);
  
  }
  real scale_to_one(real y, real start, real end) {
    real y_rescale;
     y_rescale = ((end - start) * (y - start)) / (end - start);

     return(y_rescale);
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

  int censored[S];
  real end_times;
  
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
  vector[S] mu;
  vector[S] v;
  vector[S] w;
 
  // log absolute derivative of theta_raw .* d
  for(ss in 1:S) { // jacobian adjustment is just a constant
    target += (log(fabs(fad[ss])));
    target += (log(fabs(lad[ss])));
  }
  

  // priors in general
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


  // regression statement
  for(i in 1:S) {
    scale[i] = exp(-(beta[cohort[i], 1] + 
          beta[cohort[i], 2] * occupy[i] + 
          beta[cohort[i], 3] * env[i] + 
          beta[cohort[i], 4] * (env[i]^2) +
          beta[cohort[i], 5] * leng[i])/ shape);
  }
  
  // the prior on start and end is the lifetime distribution
  dur = start - end;  // numbers head towards 0
  // this is the meat of the whole operation
  target += weibull_lpdf(dur | shape, scale);



  // assign cohort based on estimated start date
  for(i in 1:S) {
    for(j in 1:O) {
      if((start[i] <= cohort_start[j]) && 
          (start[i] > cohort_end[j])) {
        cohort[i] = j;
      }
    }
  }
 
  // conversions to make beta distribution
  for(i in 1:S) {
    mu[i] = (start[i] + mid[i] + end[i]) / 2;
    v[i] = ((mu[i] - start[i]) * (2 * mid[i] - start[i] - end[i])) / 
      ((mid[i] - mu[i]) * (end[i] - start[i]));
    w[i] = (v[i] * (end[i] - mu[i])) / (mu[i] - start[i]);
  }

  // "likelihood" of each observation, based on what species it is
  for(n in 1:N) {
    y[n] ~ pert(start[sp[n]], end[sp[n]], mid[sp[n]], 4);
  }
  
  // for those that are censored, need that last observation saying so
  for(s in 1:S) {
    if(censored[s] == 1) {
      {
        real big_e;
        big_e = scale_to_one(end_times, start[s], end[s]);
        big_e ~ pert(start[s], end[s], mid[s], 4);
      }
    }
  }
}
