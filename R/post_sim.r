library(rstan)
library(plyr)

RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(420)
nsim <- 1000

# misc functions
# poisson/negative binomial residuals


# empirical data sent to stan
souce('../R/mung.r')  # raw data
emp <- read_rdump('../data/data_dump/count_info.data.R')

# csv of marginal posteriors from stan
np <- list.files('../data/mcmc_out/', pattern = 'samp', full.names = TRUE)
preserve.fit <- read_stan_csv(np)


# for ii in seq(nsim))
#   extract random draw from marginal posteriors
#   generate fake data given known covariates
