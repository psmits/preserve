library(plyr)
library(rstan)
# processing coda into stan-like posterior draws
process.coda <- function(post, data) {
  post <- Reduce(rbind, post)
  vars <- unique(llply(str_split(colnames(post), pattern = '\\['), 
                       function(x) x[1]))

  # now break it apart
  new.post <- list()
  for(ii in seq(length(vars))) {
    temp <- post[, str_detect(colnames(post), 
                              paste0(vars[[ii]], '\\['))]
    dim.check <- str_count(colnames(temp), pattern = '\\,')
    col.num <- ncol(temp) / data$nprov
   
    if(sum(temp) == 0) {
      temp <- post[, vars[[ii]]]
      col.num <- 1
    }
    
    if(col.num == 1) {   # don't vary by time, but by region
      new.post[[ii]] <- temp
    } else if (all(dim.check <= 1)){  # vary by time, but only one variable
      holder <- array(dim = c(nrow(temp), col.num, data$nprov))
      
      colseq <- list()
      colseq[[1]] <- seq(col.num)
      for(jj in seq(from = 2, to = data$nprov)) {
        val <- colseq[[jj - 1]][length(colseq[[jj - 1]])]
        colseq[[jj]] <- seq(from = val + 1, to = val + col.num)
      }
     
      for(jj in seq(data$nprov)) {
        holder[, , jj] <- temp[, colseq[[jj]]]
      }
      new.post[[ii]] <- holder
    } else if (all(dim.check > 1)) {  # generated 3-way array z (list by region)
      holder <- list()
      regionholder <- array(dim = c(nrow(temp), col.num, data$nprov))
      
      colseq <- list()
      colseq[[1]] <- seq(col.num)
      for(jj in seq(from = 2, to = data$nprov)) {
        val <- colseq[[jj - 1]][length(colseq[[jj - 1]])]
        colseq[[jj]] <- seq(from = val + 1, to = val + col.num)
      }
     
      for(jj in seq(data$nprov)) {
        regionholder <- temp[, colseq[[jj]]]

        timeseq <- matrix(ncol = 2, nrow = data$nyear)
        timeseq[1, ] <- c(1, data$nindiv)
        for(kk in seq(from = 2, to = data$nyear))
          timeseq[kk, ] <- c((data$inidiv * (kk - 1)) + 1, data$nindiv * kk)

        timeholder <- array(dim = c(nrow(temp), data$nindiv, data$nyear))
        for(cc in seq(from = 1, to = data$nyear)) {
          timeholder[, , cc] <- regionholder[, timeseq[cc, 1]:timeseq[cc, 2]]
        }
        holder[[jj]] <- timeholder
      }
      new.post[[ii]] <- holder
    }
  }
  names(new.post) <- unlist(vars)
  new.post
}


# data for posterior predictive checks given the observed sample sizes
# post is the extraction from stan object
# data is the data used to fit the initial model
#   P: number of provinces 
#   prov: vector of province "start points"
#   C: number of time points
#   sight: sighting record (just need initial, i think)
posterior.turnover <- function(post, data) {
  # simulate for each province
  region <- list()
  for(jj in seq(data$nprov)) {
    rand <- sample(nrow(post[[1]]), 1)
    ## diversification process
    time.step <- data$nyear

    # simulate for each taxon
    life.time <- matrix(nrow = data$nindiv, ncol = data$nyear)
    life.time[, 1] <- data$y[, 1, jj]
    for(nn in seq(data$nindiv)) {
      for(ii in 2:time.step) {
        if(life.time[nn, ii - 1] == 1) {
          life.time[nn, ii] <- rbinom(1, size = 1, 
                                      prob = post$phi[rand, ii - 1, jj])
        } else {
          life.time[nn, ii] <- rbinom(1, size = 1, 
                                      prob = post$gamma[rand, ii - 1, jj])
        }
      }
    }

    sample.time <- life.time
    # observation process
    for(nn in seq(data$nindiv)) {
      for(ii in seq(from = 2, to = time.step)) {
        sample.time[nn, ii] <- life.time[nn, ii] * 
        rbinom(1, size = 1, prob = post$p[rand, ii, jj])
      }
    }

    region[[jj]] <- list(true = life.time, obs = sample.time)
  }
  region 
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
