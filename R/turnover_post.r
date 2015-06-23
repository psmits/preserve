library(plyr)
library(rstan)

data <- read_rdump('../data/data_dump/sight_info.data.R')

pat <- 'state_space_[0-9].csv'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
fit <- read_stan_csv(outs)

# need to figure out some kind of posterior predictive check
# data generating function
#   generate birth-death sequence
#   determine which are observed or not


