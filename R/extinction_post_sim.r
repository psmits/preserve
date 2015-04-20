library(rstan)
library(plyr)
RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420
nsim <- 1000

data <- read_rdump('../data/data_dump/fauna_info.data.R')

# this is for the weibull model
pat <- 'faun_expo_[0-9].csv'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
efit <- read_stan_csv(outs)

# this is for the weibull model
pat <- 'faun_weib_[0-9].csv'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
wfit <- read_stan_csv(outs)

# data setup
coh <- c(data$cohort_unc, data$cohort_cen)
gro <- c(data$group_unc, data$group_cen)
rage <- c(data$occupy_unc, data$occupy_cen)
#envs <- c(data$env_unc, data$env_cen)
size <- c(data$size_unc, data$size_cen)
duration <- c(data$dur_unc, data$dur_cen)

# weibull model
# extract values and do posterior predictive simulations
wei.fit <- extract(wfit, permuted = TRUE)
betas <- wei.fit$beta
envs <- wei.fit$env
wr <- list()
for(ii in 1:nsim) {
  n <- data$N
  alp <- sample(wei.fit$alpha, 1)
  int <- betas[sample(nrow(betas), 1), , 1]
  bet.1 <- betas[sample(nrow(betas), 1), , 2]
  bet.2 <- betas[sample(nrow(betas), 1), , 3]
  bet.3 <- betas[sample(nrow(betas), 1), , 4]
  ee <- envs[sample(nrow(envs), 1), ]

  oo <- c()
  for(jj in seq(n)) {
    reg <- int[coh[jj]] + bet.1[coh[jj]] * rage[jj] + 
    bet.2[coh[jj]] * ee[jj] + 
    bet.3[coh[jj]] * size[jj]
    oo[jj] <- rweibull(1, shape = alp, scale = exp(-(reg) / alp))
  }

  wr[[ii]] <- oo
}


# exponential model
# extract values and do posterior predictive simulations
exp.fit <- extract(efit, permuted = TRUE)
betas <- exp.fit$beta
envs <- exp.fit$env
er <- list()
for(ii in 1:nsim) {
  n <- data$N
  int <- betas[sample(nrow(betas), 1), , 1]
  bet.1 <- betas[sample(nrow(betas), 1), , 2]
  bet.2 <- betas[sample(nrow(betas), 1), , 3]
  bet.3 <- betas[sample(nrow(betas), 1), , 4]
  ee <- envs[sample(nrow(envs), 1), ]

  oo <- c()
  for(jj in seq(n)) {
    reg <- int[coh[jj]] + bet.1[coh[jj]] * rage[jj] + 
    bet.2[coh[jj]] * ee[jj] + 
    bet.3[coh[jj]] * size[jj]
    oo[jj] <- rexp(1, rate = exp(reg))
  }

  er[[ii]] <- oo
}
