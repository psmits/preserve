library(rstan)
library(plyr)

np <- list.files('../data/mcmc_out/', pattern = 'samp', full.names = TRUE)
preserve.fit <- read_stan_csv(np)
