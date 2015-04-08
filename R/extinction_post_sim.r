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
txts <- summary(fit)[[1]]
head(txts)
all(txts[, ncol(txts)] < 1.1)

# extract values and do posterior predictive simulations
extract.fit <- extract(fit, permuted = TRUE)

coh <- c(data$cohort_unc, data$cohort_cen)
gro <- c(data$group_unc, data$group_cen)

environ <- cbind(extract.fit$x_unc, extract.fit$x_cen)
taxon <- extract.fit$group
ph <- list()
for(ii in 1:nsim) {
  n <- data$N
  inter <- sample(extract.fit$intercept, 1)
  alpha <- sample(extract.fit$alpha, 1)
  slope <- sample(extract.fit$slope, 1)
  x <- environ[sample(nrow(environ), 1), ]
  gg <- taxon[sample(nrow(taxon), 1), ]

  oo <- c()
  for(jj in seq(n)) {
    reg <- inter + gg[gro[jj]] + slope * x[jj]
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
