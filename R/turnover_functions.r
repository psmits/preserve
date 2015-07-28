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
  for(jj in seq(data$J)) {
    taxa <- list()
    # simulate for each province, for the number of taxa observed?
    rand <- sample(nrow(post[[1]]), 1)
    ## diversification process

    time.step <- data$C
    life.time <- c()
    for(nn in seq(data$R)) {
      life.time[1] <- data$sight[nn, 1, jj]
      for(ii in 2:time.step) {
        if(life.time[ii - 1] == 1) {
          life.time[[ii]] <- rbinom(1, size = 1, 
                                    prob = post$phi[rand, ii - 1, jj])
        } else {
          life.time[ii] <- rbinom(1, size = 1, 
                                  prob = post$gamma[rand, ii - 1, jj])
        }
      }
    }

    sample.time <- life.time
    # observation process
    for(ii in seq(from = 2, to = time.step)) {
      sample.time[ii] <- life.time[ii] * rbinom(1, size = 1, 
                                                prob = post$p[rand, ii, jj])
    }

    taxa[[jj]] <- list(true = life.time, observed = sample.time)
  }
  true <- Reduce(rbind, llply(taxa, function(x) x[[1]]))
  observed <- Reduce(rbind, llply(taxa, function(x) x[[2]]))
  out[[kk]] <- list(true = true, observed = observed)
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
