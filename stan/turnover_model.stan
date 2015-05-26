functions {
  int first_capture(int[] y_i) {
    for (k in 1:size(y_i))
      if (y_i[k])
        return k;
    return 0;
  }
  int last_capture(int[] y_i) {
    for (k_rev in 0:(size(y_i) - 1)) {
      int k;
      k <- size(y_i) - k_rev;
      if (y_i[k])
        return k;
    }
    return 0;
  }
  vector prob_uncaptured(int T, vector p, vector phi) {
    vector[T] chi;
    chi[T] <- 1.0;
    for (t in 1:(T - 1)) {
      int t_curr;
      int t_next;
      t_curr <- T - t;
      t_next <- t_curr + 1;
      chi[t_curr] <- (1 - phi[t_curr]) + phi[t_curr] * 
      (1 - p[t_next]) * chi[t_next];
    }
    return chi;
  }
}
data {
  int N;
  int T;
  int sight[N, T];
}
transformed data {
  int<lower=0,upper=T> first[N];
  int<lower=0,upper=T> last[N];
  for (i in 1:N)
    first[i] <- first_capture(sight[i]);
  for (i in 1:N)
    last[i] <- last_capture(sight[i]);
}
parameters {
  vector<lower=0,upper=1>[T-1] phi;
  vector<lower=0,upper=1>[T] p;
}
transformed parameters {
  vector<lower=0,upper=1>[T] chi;
  chi <- prob_uncaptured(T, p, phi);
}
model {
  for(n in 1:N) {
    if(first[n] > 0) {
      for(t in first[n] + 1:last[n]) {
        1 ~ bernoulli(phi[t - 1]);
        sight[n, t] ~ bernoulli(p[t]);
      }
      1 ~ bernoulli(chi[last[n]]);
    }
  }
}
generated quantities {
  real beta;

  beta <- phi[T-1] * p[T];
}
