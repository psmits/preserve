library(plyr)
library(rstan)

# data for posterior predictive checks given the observed sample sizes
# post is the extraction from stan object
# data is the data used to fit the initial model
#   P: number of provinces 
#   prov: vector of province "start points"
#   C: number of time points
#   sight: sighting record (just need initial, i think)
posterior.turnover <- function(post, data) {
  out <- list()
  for(kk in seq(data$P)) {
    taxa <- list()
    if(kk == 1) {
      strt = 1
    } else {
      strt = data$prov[kk - 1] + 1
    }
    for(jj in seq(from = strt, to = data$prov[kk])) {
      # simulate for each province, for the number of taxa observed?
      rand <- sample(length(post$lp__), 1)
      ## diversification process
      time.step <- data$C
      life.time <- c()
      life.time[1] <- data$sight[jj, 1]
      for(ii in 2:time.step) {
        if(life.time[ii - 1] == 1) {
          life.time[[ii]] <- rbinom(1, size = 1, 
                                    prob = post$phi[rand, kk, ii - 1])
        } else if (any(life.time == 1) & life.time[ii - 1] == 0) {
          life.time[ii] <- 0
        } else {
          life.time[ii] <- rbinom(1, size = 1, 
                                  prob = post$gamma[rand, kk, ii])
        }
      }

      sample.time <- life.time
      # observation process
      for(ii in seq(from = 2, to = time.step)) {
        sample.time[ii] <- life.time[ii] * rbinom(1, size = 1, 
                                                  prob = post$p[rand, kk, ii])
      }

      taxa[[jj]] <- list(true = life.time, observed = sample.time)
    }
    true <- Reduce(rbind, llply(taxa, function(x) x[[1]]))
    observed <- Reduce(rbind, llply(taxa, function(x) x[[2]]))
    out[[kk]] <- list(true = true, observed = observed)
  }
  out
}


# count how many taxa actually seen from the simulation
# out is just a single provincial simulation
# set can be 'true' or 'observed'
check.count <- function(out, set = 'true') {
  chck <- c()
  if(set == 'true') {
    for(ii in seq(nrow(out$true))) {
      chck[ii] <- any(out$true[ii, ] != 0)
    }
  } else if (set == 'observed') {
    for(ii in seq(nrow(out$true))) {
      chck[ii] <- any(out$observed[ii, ] != 0)
    }
  }
  chck
}
