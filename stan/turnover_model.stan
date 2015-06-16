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
  real state_space_log(int[] y, vector phi, vector p, vector gamma) {
    int ft;
    int lt;
    int S;
    vector[first_capture(y) * (size(y) - last_capture(y) + 1)] lp;
    int i;

    ft <- first_capture(y);
    lt <- last_capture(y);
    S <- size(y);
    i <- 1;

    for(t_first_alive in 1:ft) {
      for (t_last_alive in lt:S) {
        real sl;
        int z[S];

        for(l in 1:S) {
          z[l] <- 0;
        }
        for(a in t_first_alive:t_last_alive) {
          z[a] <- 1;
        }

        sl <- bernoulli_log(z[1], gamma[1]);
        for(j in 2:S) {
          sl <- sl + bernoulli_log(z[j], (z[j - 1] * phi[j - 1]) + 
              (gamma[j] * (1 - z[j - 1])));
        }
        for(k in 1:S) {
          sl <- sl + bernoulli_log(y[k], z[k] * p[k]);
        }

        lp[i] <- sl;
        i <- i + 1;
      }
    }
    return log_sum_exp(lp);
  }
}
data {
  int R;  // rows, occurrences records
  int C;  // columns, temporal bins
  int T;  // number of unique taxa
  int P;  // number of provinces
  int sight[R, C];  // total sight record matrix
  int taxon[R];  // which taxon for record
  int prov[P];  // size of each province
}
parameters {
  vector<lower=0,upper=1>[C - 1] phi[P];
  vector<lower=0,upper=1>[C] p[P];
  vector<lower=0,upper=1>[C] gamma[P];
}
model {
  for(k in 1:P) {
    if(k == 1) {
      for(y in 1:(prov[k])) {
        sight[y] ~ state_space(phi[k], p[k], gamma[k]);
      }
    } else {
      for(y in (prov[k - 1] + 1):prov[k]) {
        sight[y] ~ state_space(phi[k], p[k], gamma[k]);
      }
    }
  }
}
