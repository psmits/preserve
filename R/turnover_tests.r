library(rstan)
library(plyr)
RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420
nsim <- 100

# data generating function
generate.capture <- function(tt = 30, surv = 0.3, gg = 0.3, samp = 0.3) {
  present <- c()
  run.prob <- c()
  alive <- c()
  k <- 1
  for(ii in seq(tt)) {
    if(ii == 1) {
      run.prob[ii] <- gg
    } else if(ii > 1) {
      k <- k * (1 - run.prob[ii - 1])
      run.prob[ii] <- run.prob[ii - 1] * surv + gg * k
    }
    alive[ii] <- runif(1) <= run.prob[ii] * samp
  }
  out <- list(alive = as.numeric(alive), 
              run.prob = run.prob)
  out
}

sims <- replicate(nsim, expr = {generate.capture()}, simplify = FALSE)
record <- Reduce(rbind, llply(sims, function(x) x$alive))

num <- nrow(record)
tim <- ncol(record)
data <- list(N = num, T = tim, sight = record)

with(data, {stan_rdump(list = c('N', 'T', 'sight'),
                       file = '../data/data_dump/sim_info.data.R')})
# then fit model

# after model fit
#pat <- 'state_space_[0-9].csv'
#outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
#fit <- read_stan_csv(outs)
