library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(arm)
library(rstan)
library(survival)
library(stringr)
library(grid)
library(gridBase)
library(gridExtra)
library(xtable)
library(ellipse)
library(loo)
library(ppcor)
source('../R/mung.r')
source('../R/multiplot.r')
source('../R/borrow_plotcorr.r')
set.seed(420)

data.file <- list.files('../data', pattern = 'Occs')
fossil <- read.csv(paste0('../data/', data.file))
bibr <- fossil
payne <- read.table('../data/payne_bodysize/Occurrence_PaleoDB.txt',
                    header = TRUE, stringsAsFactors = FALSE)

lump.file <- list.files('../data', pattern = 'lump')
lump <- read.csv(paste0('../data/', lump.file))
gts <- rev(as.character(lump[, 2]))

sepkoski.data <- sort.data(bibr, payne, taxon = 'Rhynchonellata', 
                           bins = 'StageNewOrdSplitNoriRhae20Nov2013', 
                           gts = gts,
                           cuts = 'Chang',
                           bot = 'Trem')

data <- read_rdump('../data/data_dump/impute_info.data.R')

pat <- 'faun_'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
ids <- rep(1:(length(outs) / 4), each = 4)
outs <- split(outs, ids)


wfit <- read_stan_csv(outs[[2]])
log.lik <- extract_log_lik(wfit)
wei.fit <- rstan::extract(wfit, permuted = TRUE)
# this will need to be updated with number of models
npred <- 7


# environmental preference
prob.epi <- qbeta(sepkoski.data$epi / (sepkoski.data$epi + sepkoski.data$off), 
                  sepkoski.data$epi.bck + 1, sepkoski.data$off.bck + 1)
env.odds <- rescale(prob.epi)
# start asking questions and making tables


## partial correlation between duration, occupancy, preservation
#par.cor <- data.frame(dur = sepkoski.data$duration, 
#                      oc = rescale(logit(sepkoski.data$occupy)), 
#                      gp = rescale(sepkoski.data$gap))
#cor.mat <- cor(par.cor)            # correlation matrix
#pcor.mat <- pcor(par.cor)$estimate  # partial correlation matrix
#
#n <- c('duration', 'geographic range', 'gap ratio')
#rownames(cor.mat) <- colnames(cor.mat) <- n
#rownames(pcor.mat) <- colnames(pcor.mat) <- n
#diag(cor.mat) <- diag(pcor.mat) <- NA
#cor.tex <- xtable(cor.mat, label = 'tab:corr')
#pcor.tex <- xtable(pcor.mat, label = 'tab:pcor')
#print.xtable(cor.tex, file = '../doc/covariate_cor.tex')
#print.xtable(pcor.tex, file = '../doc/covariate_pcor.tex')

oo <- log(sepkoski.data$occupy)
sd(oo)
(max(oo) - min(oo)) / sd(oo)

sd(prob.epi)

td <- apply(wei.fit$samp, 1, sd)

dd <- wei.fit$delta * 2 * td

sum(abs(wei.fit$mu_prior[, 2]) > dd) / length(dd)
sum(abs(wei.fit$mu_prior[, 3]) > dd) / length(dd)
sum(abs(wei.fit$mu_prior[, 4]) > dd) / length(dd)


tt <- (abs((wei.fit$mu_prior[, 2] * 2 * sd(oo))) - wei.fit$delta) > 0
et <- (abs((wei.fit$mu_prior[, 3] * 2 * sd(prob.epi))) - wei.fit$delta) < 0
sum(tt) / length(tt)
sum(et) / length(et)

et <- (abs((wei.fit$beta[, 17, 3] * 2 * sd(prob.epi))) - wei.fit$delta) > 0
sum(et) / length(et)
# cohort 16 17 33



# mean of all coefficients
# sd of all coefficients
name.mu <- c('mu_i', 'mu_r', 'mu_v', 'mu_v2', 'mu_rxv', 'mu_rxv2', 'mu_m', 'mu_s')
name.tau <- c('tau_i', 'tau_r', 'tau_v', 'tau_v2', 'tau_rxv', 
              'tau_rxv2', 'tau_m', 'tau_s')
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
print.xtable(param.tex, file = '../doc/table_param.tex')


# tail probabilities
p.r <- 1 - (pnorm(0, wei.fit$mu_prior[, 2], wei.fit$sigma[, 2]))
mean(p.r)
sd(p.r)

p.v <- 1 - (pnorm(0, wei.fit$mu_prior[, 3], wei.fit$sigma[, 3]))
mean(p.v)
sd(p.v)

p.v2 <- (pnorm(0, wei.fit$mu_prior[, 4], wei.fit$sigma[, 4]))
mean(p.v2)
sd(p.v2)


# inflection point is defined - (beta_v) / 2(beta_v2)
# probability that average of the parabola defined
#   log(sigma) = -(mu_i + mu_v * v + mu_v2 * v^2) / alpha
# midpoint favors epicontinental
#coef.inf <- wei.fit$mu_prior[, 3:4]
#inf.rel <- sum((-(coef.inf[, 1])) / (2 * coef.inf[, 2]) > 0) / 4000
#
#inf.cohort <- inf.low <- inf.hgh <- c()
#
#for(ii in seq(data$O)) {
#  h <- wei.fit$beta[, ii, 3:4]
#  inf.cohort[ii] <- sum((-(h[, 1])) / (2 * h[, 2]) > 0) / 4000
#  inf.low[ii] <- sum((-(h[, 1])) / (2 * h[, 2]) > 0.5) / 4000
#  inf.hgh[ii] <- sum((-(h[, 1])) / (2 * h[, 2]) < (-0.5)) / 4000
#}
#inf.cohort[6:15]
#
#
## probabilty of negative correlation term
#cor.int.range <- sum((wei.fit$Omega[, 1, 2] < 0)) / 4000
#cor.int.env <- sum((wei.fit$Omega[, 1, 3] < 0)) / 4000
#cor.range.env <- sum((wei.fit$Omega[, 2, 3] > 0)) / 4000
#
#
#
## but but but
##   massive difference in the between cohort variances of these two regression coefficients
##   see param.est table!
#tv.gt.sr <- sum(apply(wei.fit$sigma[, 2:3], 1, function(x) x[2] > x[1])) / 4000
#
