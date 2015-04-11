library(rstan)
library(plyr)
RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420
nsim <- 1000

data <- read_rdump('../data/data_dump/fauna_info.data.R')

# this is for the total model
pat <- paste0('faun_surv_', '[0-9].csv')
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
fit <- read_stan_csv(outs)
txt <- summary(fit)[[1]]
all(txt[, ncol(txt)] < 1.1, na.rm = TRUE)

# extract values and do posterior predictive simulations
extract.fit <- extract(fit, permuted = TRUE)

coh <- c(data$cohort_unc, data$cohort_cen)
gro <- c(data$group_unc, data$group_cen)
rage <- c(data$occupy_unc, data$occupy_cen)
envs <- c(data$env_unc, data$env_cen)

environ <- cbind(extract.fit$x_unc, extract.fit$x_cen)
coefs <- extract.fit$coef
ph <- list()
for(ii in 1:nsim) {
  n <- data$N
  alpha <- sample(extract.fit$alpha, 1)
  int <- coefs[sample(nrow(coefs), 1), , 1]
  slp <- coefs[sample(nrow(coefs), 1), , 2]

  oo <- c()
  for(jj in seq(n)) {
    reg <- int[coh[jj]] + slp[coh[jj]] * rage[jj]
    oo[jj] <- rweibull(1, scale = exp(-(reg) / alpha), shape = alpha)
  }

  ph[[ii]] <- oo
}

duration <- c(data$dur_unc, data$dur_cen)
#par(mfrow = c(5, 4), mar = c(4, 4, 2, 2))
#hist(duration, xlab = '', main = 'y')
#for(s in 1:19) hist(ph[[s]], xlab = '', main = paste('y_rep', s))
tstat.mean <- sum(laply(ph, mean) > mean(duration))
tstat.med <- sum(laply(ph, median) > median(duration))
tstat.75 <- sum(laply(ph, function(x) quantile(x, .75)) > quantile(duration, .75))
tstat.25 <- sum(laply(ph, function(x) quantile(x, .25)) > quantile(duration, .25))
