library(rstan)
library(plyr)

RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(420)
nsim <- 1000

# misc functions
# poisson/negative binomial residuals


# empirical data sent to stan
source('../R/mung.r')  # raw data
emp <- read_rdump('../data/data_dump/count_info.data.R')

# csv of marginal posteriors from stan
np <- list.files('../data/mcmc_out/', pattern = 'samp', full.names = TRUE)
preserve.fit <- read_stan_csv(np)

for(ii in seq(nsim)) {
  n <- emp$C
  off <- emp$off
  
  phi.est <- sample(preserve.fit$phi, 1)

  # for each observation estimate the sighting
  # get genus membership for each observation
  # get that genus posterior estimate

  oo <- c()
  for(jj in seq(n)) {
    #oo[jj] <- rnbinom(n = 1, mu = , size = phi.est)
  }
}
