library(rstan)
library(plyr)
RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420
nsim <- 1000

#map <- TRUE
#error <- FALSE

# functions for calculating the residuals
wei.surv <- function(time, scale, shape) {
  exp(-((1 / scale) * time)**shape)
}

cum.haz <- function(time, scale, shape) {
  -log(wei.surv(time, scale, shape))
}

martingale.res <- function(time, scale, shape, inclusion) {
  inclusion - cum.haz(time, scale, shape)
}

deviance.res <- function(time, scale, shape, inclusion) {
  martin <- martingale.res(time, scale, shape, inclusion)
  si <- sign(martin)
  inn <- martin + (inclusion * log(inclusion - martin))
  under <- -2 * inn
  out <- si * sqrt(under)
  out
}

# bring in data
data <- read_rdump('../data/data_dump/fauna_info.data.R')

if(map) {
  # for MAP models
  # this is for the weibull model
  pat <- 'faun_expo_map_[0-9].csv'
  outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
  efit <- read_stan_csv(outs)

  #a <- summary(efit)[[1]]
  #all(a[, ncol(a)] < 1.1, na.rm = TRUE)

  # this is for the weibull model
  pat <- 'faun_weib_map_[0-9].csv'
  outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
  wfit <- read_stan_csv(outs)

  #a <- summary(wfit)[[1]]
  #all(a[, ncol(a)] < 1.1, na.rm = TRUE)

  # data setup
  coh <- c(data$cohort_unc, data$cohort_cen)
  gro <- c(data$group_unc, data$group_cen)
  rage <- c(data$occupy_unc, data$occupy_cen)
  envs <- c(data$env_unc, data$env_cen)
  lits <- c(data$lit_unc, data$lit_cen)
  size <- c(data$size_unc, data$size_cen)
  duration <- c(data$dur_unc, data$dur_cen)
  extinct <- c(rep(1, data$N_unc), rep(0, data$N_cen))

  # weibull model
  # extract values and do posterior predictive simulations
  wei.fit <- extract(wfit, permuted = TRUE)
  betas <- wei.fit$beta
  wr <- list()
  wr.res <- list()
  for(ii in 1:nsim) {
    n <- data$N
    alp <- sample(wei.fit$alpha, 1)
    int <- betas[sample(nrow(betas), 1), , 1]
    bet.1 <- betas[sample(nrow(betas), 1), , 2]
    bet.2 <- betas[sample(nrow(betas), 1), , 3]
    bet.3 <- betas[sample(nrow(betas), 1), , 4]
    bet.4 <- betas[sample(nrow(betas), 1), , 5]

    oo <- c()
    for(jj in seq(n)) {
      reg <- int[coh[jj]] + bet.1[coh[jj]] * rage[jj] + 
      bet.2[coh[jj]] * envs[jj] + 
      bet.3[coh[jj]] * (envs[jj] * envs[jj]) + 
      bet.4[coh[jj]] * size[jj]
      oo[jj] <- rweibull(1, shape = alp, scale = exp(-(reg) / alp))
      rr[jj] <- deviance.res(duration[jj], 
                            scale = exp(reg), 
                            shape = alp, 
                            inclusion = extinct[jj])
    }

    wr[[ii]] <- oo
    wr.res[[ii]] <- rr
  }

  # exponential model
  # extract values and do posterior predictive simulations
  exp.fit <- extract(efit, permuted = TRUE)
  betas <- exp.fit$beta
  er <- list()
  er.res <- list()
  for(ii in 1:nsim) {
    n <- data$N
    int <- betas[sample(nrow(betas), 1), , 1]
    bet.1 <- betas[sample(nrow(betas), 1), , 2]
    bet.2 <- betas[sample(nrow(betas), 1), , 3]
    bet.3 <- betas[sample(nrow(betas), 1), , 4]
    bet.4 <- betas[sample(nrow(betas), 1), , 5]

    oo <- c()
    for(jj in seq(n)) {
      reg <- int[coh[jj]] + bet.1[coh[jj]] * rage[jj] + 
      bet.2[coh[jj]] * envs[jj] + 
      bet.3[coh[jj]] * (envs[jj] * envs[jj]) + 
      bet.4[coh[jj]] * size[jj]
      oo[jj] <- rexp(1, rate = exp(reg))
      rr[jj] <- deviance.res(duration[jj], 
                            scale = 1 / exp(reg), 
                            shape = 1, 
                            inclusion = extinct[jj])
    }

    er[[ii]] <- oo
    er.res[[ii]] <- rr
  }
} else {
  # for measurement error models
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
  wr.res <- list()
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
      rr[jj] <- deviance.res(duration[jj], 
                            scale = exp(reg), 
                            shape = alp, 
                            inclusion = extinct[jj])
    }

    wr[[ii]] <- oo
    wr.res[[ii]] <- rr
  }


  # exponential model
  # extract values and do posterior predictive simulations
  exp.fit <- extract(efit, permuted = TRUE)
  betas <- exp.fit$beta
  envs <- exp.fit$env
  er <- list()
  er.res <- list()
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
      rr[jj] <- deviance.res(duration[jj], 
                            scale = 1 / exp(reg), 
                            shape = 1, 
                            inclusion = extinct[jj])
    }

    er[[ii]] <- oo
    er.res[[ii]] <- rr
  }
}
