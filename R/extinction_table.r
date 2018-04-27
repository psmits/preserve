library(tidyverse)
library(arm)
library(rstan)
library(xtable)

data <- read_rdump('../data/data_dump/impute_info.data.R')

pat <- 'surv_cweib_cens'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)

wfit <- read_stan_csv(outs)
#log.lik <- extract_log_lik(wfit)
wei.fit <- rstan::extract(wfit, permuted = TRUE)
# this will need to be updated with number of models
npred <- 5


# mean of all coefficients
# sd of all coefficients
name.mu <- c('mu_i', 'mu_r', 'mu_v', 'mu_v2', 'mu_m', 'mu_s')
name.tau <- c('tau_i', 'tau_r', 'tau_v', 'tau_v2', 'tau_m')
qq <- c(0.1, 0.5, 0.9)
param.est <- rbind(data.frame(p = name.mu[1:npred],
                              m = apply(wei.fit$mu_prior, 2, mean), 
                              s = apply(wei.fit$mu_prior, 2, sd),
                              qa = t(apply(wei.fit$mu_prior, 2, 
                                           function(x) quantile(x, qq)))), 
                   data.frame(p = name.tau[1:npred],
                              m = apply(wei.fit$sigma, 2, mean), 
                              s = apply(wei.fit$sigma, 2, sd),
                              qa = t(apply(wei.fit$sigma, 2, 
                                     function(x) quantile(x, qq)))),
                   data.frame(p = 'delta',
                              m = mean(wei.fit$delta),
                              s = sd(wei.fit$delta),
                              qa = t(quantile(wei.fit$delta, qq))),
                   data.frame(p = 'alpha',
                              m = mean(wei.fit$alpha),
                              s = sd(wei.fit$alpha),
                              qa = t(quantile(wei.fit$alpha, qq))))

param.tex <- xtable(param.est, label = 'tab:param')
print.xtable(param.tex, file = '../doc/table_param_draft.tex')


# tail probabilities
p.r <- (pnorm(0, wei.fit$mu_prior[, 2], wei.fit$sigma[, 2]))

p.v <- 1 - (pnorm(0, wei.fit$mu_prior[, 3], wei.fit$sigma[, 3]))
mean(p.v)
sd(p.v)


p.v2 <- (pnorm(0, wei.fit$mu_prior[, 4], wei.fit$sigma[, 4]))
mean(p.v2)
sd(p.v2)

# inflection point is defined - (beta_v) / 2(beta_v2)
# probability that average of the parabola defined
#   log(sigma) = -(mu_i + mu_v * v + mu_v2 * v^2) / alpha
coef.inf <- wei.fit$mu_prior[, 3:4]
inf.rel <- sum((-(coef.inf[, 1])) / (2 * coef.inf[, 2]) > 0) / 4000
# midpoint favors epicontinental

inf.cohort <- inf.low <- inf.hgh <- c()
for(ii in seq(data$O)) {
  h <- wei.fit$beta[, ii, 3:4]
  inf.cohort[ii] <- sum((-(h[, 1])) / (2 * h[, 2]) > 0) / 4000
  inf.low[ii] <- sum((-(h[, 1])) / (2 * h[, 2]) > 0.1) / 4000
  inf.hgh[ii] <- sum((-(h[, 1])) / (2 * h[, 2]) < (-0.1)) / 4000
}
#inf.cohort # cohort pref epi
#inf.low # strongish for epi??
#inf.hgh # strongish against epi?


# probabilty of correlations in some directions
cor.int.range <- sum((wei.fit$Omega[, 1, 2] < 0)) / 4000 # p(negative)
cor.int.env <- sum((wei.fit$Omega[, 1, 3] < 0)) / 4000 # p(negative)
cor.int.env2 <- sum((wei.fit$Omega[, 1, 4] < 0)) / 4000 # p(negative)
cor.range.env <- sum((wei.fit$Omega[, 2, 3] > 0)) / 4000 # p(positive)
cor.range.env2 <- sum((wei.fit$Omega[, 2, 4] > 0)) / 4000 # p(positive)


# but but but
#   massive difference in the between cohort variances of these two regression coefficients
#   see param.est table!
tv.gt.sr <- sum(apply(wei.fit$sigma[, 2:3], 1, function(x) x[2] > x[1])) / 4000
