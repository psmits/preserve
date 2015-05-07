library(rstan)
library(plyr)
RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420
nsim <- 1000

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

#pat <- 'faun_expo_map_[0-9].csv'
#outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
#efit <- read_stan_csv(outs)
#pat <- 'faun_weib_map_[0-9].csv'
#outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
#wfit <- read_stan_csv(outs)
#pat <- 'faun_expo_[0-9].csv'
#outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
#efit <- read_stan_csv(outs)
#wfit <- read_stan_csv(outs)

post.sim <- function(data, fit, map = FALSE, expo = TRUE) {
  # data setup
  coh <- c(data$cohort_unc, data$cohort_cen)
  gro <- c(data$group_unc, data$group_cen)
  rage <- c(data$occupy_unc, data$occupy_cen)
  envs <- c(data$env_unc, data$env_cen)
  lits <- c(data$lit_unc, data$lit_cen)
  size <- c(data$size_unc, data$size_cen)
  duration <- c(data$dur_unc, data$dur_cen)
  extinct <- c(rep(1, data$N_unc), rep(0, data$N_cen))

  if(map) {
    if(expo) {
      exp.fit <- rstan::extract(fit, permuted = TRUE)
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
        rr <- c()
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
      out <- list(er, er.res)
    } else if(!expo) {
      # weibull model
      # extract values and do posterior predictive simulations
      wei.fit <- rstan::extract(fit, permuted = TRUE)
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
        rr <- c()
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
      out <- list(wr, wr.res)
    }
  } else if(!map) {
    # for measurement error models
    if(expo) {
      exp.fit <- rstan::extract(fit, permuted = TRUE)
      betas <- exp.fit$beta
      envs <- exp.fit$env
#      lits <- exp.fit$lit
      er <- list()
      er.res <- list()
      for(ii in 1:nsim) {
        n <- data$N
        int <- betas[sample(nrow(betas), 1), , 1]
        bet.1 <- betas[sample(nrow(betas), 1), , 2]
        bet.2 <- betas[sample(nrow(betas), 1), , 3]
        bet.3 <- betas[sample(nrow(betas), 1), , 4]
        bet.4 <- betas[sample(nrow(betas), 1), , 5]
#        bet.5 <- betas[sample(nrow(betas), 1), , 6]
#        bet.6 <- betas[sample(nrow(betas), 1), , 7]
#        bet.7 <- betas[sample(nrow(betas), 1), , 8]
#        bet.8 <- betas[sample(nrow(betas), 1), , 9]
        ee <- envs[sample(nrow(envs), 1), ]
#        ll <- lits[sample(nrow(lits), 1), ]

        oo <- c()
        rr <- c()
        for(jj in seq(n)) {
          reg <- int[coh[jj]] + bet.1[coh[jj]] * rage[jj] + 
          bet.2[coh[jj]] * ee[jj] + 
          bet.3[coh[jj]] * (ee[jj])^2
          bet.4[coh[jj]] * size[jj]
#          + bet.5[coh[jj]] * ll[jj] + 
#          bet.6[coh[jj]] * (ll[jj])^2 + 
#          bet.7[coh[jj]] * (ll[jj] * ee[jj]) + 
#          bet.8[coh[jj]] * (ll[jj] * ee[jj])^2

          oo[jj] <- rexp(1, rate = exp(reg))
          rr[jj] <- deviance.res(duration[jj], 
                                 scale = 1 / exp(reg), 
                                 shape = 1, 
                                 inclusion = extinct[jj])
        }

        er[[ii]] <- oo
        er.res[[ii]] <- rr
      }
      out <- list(er, er.res)
    } else if(!expo) {

      # weibull model
      # extract values and do posterior predictive simulations
      wei.fit <- rstan::extract(fit, permuted = TRUE)
      betas <- wei.fit$beta
      envs <- wei.fit$env
#      lits <- wei.fit$lit
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
#        bet.5 <- betas[sample(nrow(betas), 1), , 6]
#        bet.6 <- betas[sample(nrow(betas), 1), , 7]
#        bet.7 <- betas[sample(nrow(betas), 1), , 8]
#        bet.8 <- betas[sample(nrow(betas), 1), , 9]
        ee <- envs[sample(nrow(envs), 1), ]
#        ll <- lits[sample(nrow(lits), 1), ]

        oo <- c()
        rr <- c()
        for(jj in seq(n)) {
          reg <- int[coh[jj]] + bet.1[coh[jj]] * rage[jj] + 
          bet.2[coh[jj]] * ee[jj] + 
          bet.3[coh[jj]] * (ee[jj])^2
          bet.4[coh[jj]] * size[jj]
#          + bet.5[coh[jj]] * ll[jj] + 
#          bet.6[coh[jj]] * (ll[jj])^2 + 
#          bet.7[coh[jj]] * (ll[jj] * ee[jj]) + 
#          bet.8[coh[jj]] * (ll[jj] * ee[jj])^2
          oo[jj] <- rweibull(1, shape = alp, scale = exp(-(reg) / alp))
          rr[jj] <- deviance.res(duration[jj], 
                                 scale = exp(reg), 
                                 shape = alp, 
                                 inclusion = extinct[jj])
        }

        wr[[ii]] <- oo
        wr.res[[ii]] <- rr
      }
      out <- list(wr, wr.res)
    }
  }
  return(out)
}

