library(plyr)
library(rstan)
source('../R/turnover_functions.r')

data <- read_rdump('../data/data_dump/sight_info.data.R')

pat <- 'state_space_[0-9].csv'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
fit <- read_stan_csv(outs)
post <- extract(fit, permuted = TRUE)
# 3-way arrays
#   [sample, location, time-point]

post.check <- replicate(100, posterior.turnover(post = post, data = data), 
                        simplify = FALSE)

true.seen <- obs.seen <- div.true <- div.obs <- list()
for(ii in seq(length(post.check))) {
  true.seen[[ii]] <- lapply(post.check[[ii]], function(x) 
                            check.count(x, set = 'true'))
  obs.seen[[ii]] <- lapply(post.check[[ii]], function(x) 
                           check.count(x, set = 'observed'))

  div.true[[ii]] <- llply(post.check[[ii]], function(x) colSums(x$true))
  div.obs[[ii]] <- llply(post.check[[ii]], function(x) colSums(x$obs))
}
div.dist <- llply(1:4, function(y) 
                  Reduce(rbind, llply(div.obs, function(x) x[[y]])))

#plot(div.dist[[1]][1, ], type = 'l')
#for(ii in 2:nrow(div.dist[[1]])) lines(div.dist[[1]][ii, ])
#lines(colSums(data$sight[1:data$prov[1], ]), col = 'blue')


# simulations with even "effort" across provinces
sim.data = list(P = 4, 
                prov = seq(from = 1000, to = 4000, by = 1000), 
                C = data$C,
                sight = matrix(rep(0, 4000), ncol = 1))
post.sim <- replicate(100, posterior.turnover(post = post, data = sim.data),
                        simplify = FALSE)
true.seen <- obs.seen <- div.true <- div.obs <- list()
for(ii in seq(length(post.sim))) {
  true.seen[[ii]] <- lapply(post.sim[[ii]], function(x) 
                            check.count(x, set = 'true'))
  obs.seen[[ii]] <- lapply(post.sim[[ii]], function(x) 
                           check.count(x, set = 'observed'))

  div.true[[ii]] <- llply(post.sim[[ii]], function(x) colSums(x$true))
  div.obs[[ii]] <- llply(post.sim[[ii]], function(x) colSums(x$obs))
}
# i want to see the "true" curves
div.sim <- llply(1:4, function(y) 
                 Reduce(rbind, llply(div.true, function(x) x[[y]])))

#plot(div.dist[[1]][1, ], type = 'l')
#for(ii in 2:nrow(div.dist[[1]])) lines(div.dist[[1]][ii, ])
