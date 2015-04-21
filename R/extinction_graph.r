library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(rstan)
library(survival)
library(stringr)
library(grid)

map <- TRUE
error <- FALSE
source('../R/extinction_post_sim.r')

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 20),
             legend.text = element_text(size = 15),
             legend.title = element_text(size = 16),
             legend.key.size = unit(2, 'cm'),
             strip.text = element_text(size = 10))

# data setup
coh <- c(data$cohort_unc, data$cohort_cen)
gro <- c(data$group_unc, data$group_cen)
rage <- c(data$occupy_unc, data$occupy_cen)
envs <- c(data$env_unc, data$env_cen)  # maximum a posteriori estimate
size <- c(data$size_unc, data$size_cen)
duration <- c(data$dur_unc, data$dur_cen)

# lets make survival curves
condition <- c(rep(1, data$N_unc), rep(0, data$N_cen))
condition[duration == 1 & condition == 1] <- 2

emp.surv <- survfit(Surv(time = duration, time2 = duration, 
                         event = condition, type = 'interval') ~ 1)
emp.surv <- data.frame(cbind(time = emp.surv$time, surv = emp.surv$surv))
emp.surv <- rbind(c(0, 1), emp.surv)

# weibull model
wei.surv <- llply(wr, function(x) survfit(Surv(x) ~ 1))
wei.surv <- llply(wei.surv, function(x) {
                  y <- data.frame(time = x$time, surv = x$surv)
                  rbind(c(0, 1), y)
                  y})
wei.surv <- Reduce(rbind, Map(function(x, y) {
                              x$group <- y
                              x}, 
                              x = wei.surv, 
                              y = seq(length(wei.surv))))
wei.surv$label <- 'Weibull'

# exponential model
exp.surv <- llply(er, function(x) survfit(Surv(x) ~ 1))
exp.surv <- llply(exp.surv, function(x) {
                  y <- data.frame(time = x$time, surv = x$surv)
                  rbind(c(0, 1), y)
                  y})
exp.surv <- Reduce(rbind, Map(function(x, y) {
                              x$group <- y
                              x}, 
                              x = exp.surv, 
                              y = seq(length(exp.surv))))
exp.surv$label <- 'Exponential'

# combine both simulation sets
sim.surv <- rbind(wei.surv, exp.surv)

# fit model
surv.plot <- ggplot(emp.surv, aes(x = time, y = surv))
surv.plot <- surv.plot + geom_line(data = sim.surv, 
                                   aes(x = time, y = surv, group = group),
                                   colour = 'grey', alpha = 0.2)
surv.plot <- surv.plot + geom_line(size = 1)
surv.plot <- surv.plot + coord_cartesian(xlim = c(-0.5, max(duration) + 2))
surv.plot <- surv.plot + facet_grid(. ~ label, labeller = label_parsed)
surv.plot <- surv.plot + labs(x = 'Duration in stages', y = 'P(T > t)')
ggsave(surv.plot, filename = '../doc/survival/figure/suvival_curves.png',
       width = 8, height = 5)

# make plot of correlation and covariance matrices
# row is sample
# dim 2 is row
# dim 3 is col
get.covcor <- function(stanfit) {
  cor.median <- matrix(, ncol = 5, nrow = 5)
  cor.mean <- matrix(, ncol = 5, nrow = 5)
  cor.10 <- matrix(, ncol = 5, nrow = 5)
  cor.90 <- matrix(, ncol = 5, nrow = 5)
  cov.median <- matrix(, ncol = 5, nrow = 5)
  cov.mean <- matrix(, ncol = 5, nrow = 5)
  cov.10 <- matrix(, ncol = 5, nrow = 5)
  cov.90 <- matrix(, ncol = 5, nrow = 5)
  for(ii in seq(5)) {
    for(jj in seq(5)) {
      cor.median[ii, jj] <- median(stanfit$Omega[, ii, jj])
      cor.mean[ii, jj] <- mean(stanfit$Omega[, ii, jj])
      cor.10[ii, jj] <- quantile(stanfit$Omega[, ii, jj], probs = .1)
      cor.90[ii, jj] <- quantile(stanfit$Omega[, ii, jj], probs = .9)
      cov.median[ii, jj] <- median(stanfit$Sigma[, ii, jj])
      cov.mean[ii, jj] <- mean(stanfit$Sigma[, ii, jj])
      cov.10[ii, jj] <- quantile(stanfit$Sigma[, ii, jj], probs = .1)
      cov.90[ii, jj] <- quantile(stanfit$Sigma[, ii, jj], probs = .9)
    }
  }

  out <- list(cor.median, cor.mean, cor.10, cor.90, 
              cov.median, cov.mean, cov.10, cov.90)
  out <- llply(out, function(x) {
               rownames(x) <- c('i', 'r', 'e', 'l', 'm')
               colnames(x) <- c('i', 'r', 'e', 'l', 'm')
               x})
  out
}
wei.covcor <- get.covcor(wei.fit)
exp.covcor <- get.covcor(exp.fit)

relab.x <- scale_x_discrete(labels = c('i' = expression(beta[intercept]), 
                                       'r' = expression(beta[range]),
                                       'e' = expression(beta[environment]), 
                                       'l' = expression(beta[lithology]), 
                                       'm' = expression(beta[size])))
relab.y <- scale_y_discrete(labels = c('i' = expression(beta[intercept]), 
                                       'r' = expression(beta[range]),
                                       'e' = expression(beta[environment]), 
                                       'l' = expression(beta[lithology]), 
                                       'm' = expression(beta[size])))
# correlation matrix
omega.med <- melt(list(Exponential = exp.covcor[[1]], 
                       Weibull = wei.covcor[[1]]))
omega.plot <- ggplot(omega.med, aes(x = Var1, y = Var2, fill = value))
omega.plot <- omega.plot + geom_tile()
omega.plot <- omega.plot + facet_grid(. ~ L1, labeller = label_parsed)
omega.plot <- omega.plot + scale_fill_gradient2(name = 'Median\nCorrelation',
                                                low = 'blue', 
                                                mid = 'white', 
                                                high = 'red')
omega.plot <- omega.plot + relab.x + relab.y
omega.plot <- omega.plot + labs(x = '', y = '')
ggsave(omega.plot, filename = '../doc/survival/figure/correlation_heatmap.png',
       width = 10, height = 5)

# covariance matrix
sigma.med <- melt(list(Exponential = exp.covcor[[5]], 
                       Weibull = wei.covcor[[5]]))
sigma.plot <- ggplot(sigma.med, aes(x = Var1, y = Var2, fill = value))
sigma.plot <- sigma.plot + geom_tile()
sigma.plot <- sigma.plot + facet_grid(. ~ L1, labeller = label_parsed)
sigma.plot <- sigma.plot + scale_fill_gradient2(name = 'Median\nCovariance', 
                                                low = 'blue', 
                                                mid = 'white', 
                                                high = 'red')
sigma.plot <- sigma.plot + relab.x + relab.y
sigma.plot <- sigma.plot + labs(x = '', y = '')
ggsave(sigma.plot, filename = '../doc/survival/figure/covariance_heatmap.png',
       width = 10, height = 5)

# histogram of posterior of correlation between inter and env
baseline.covar <- data.frame(value = c(exp.fit$Omega[, 1, 3], 
                                       wei.fit$Omega[, 1, 3],
                                       exp.fit$Omega[, 1, 4], 
                                       wei.fit$Omega[, 1, 4],
                                       exp.fit$Omega[, 1, 5], 
                                       wei.fit$Omega[, 1, 5]),
                             lab = c(rep('Exponential', 
                                         length(exp.fit$Omega[, 1, 3])),
                                     rep('Weibull', 
                                         length(exp.fit$Omega[, 1, 3])),
                                     rep('Exponential', 
                                         length(wei.fit$Omega[, 1, 3])),
                                     rep('Weibull', 
                                         length(wei.fit$Omega[, 1, 3])),
                                     rep('Exponential', 
                                         length(wei.fit$Omega[, 1, 3])),
                                     rep('Weibull', 
                                         length(wei.fit$Omega[, 1, 3]))))
baseline.covar$var <- c(rep('Cor(beta[intercept], beta[environment])',
                            2 * length(exp.fit$Omega[, 1, 3])),
                        rep('Cor(beta[intercept], beta[lithology])',
                            2 * length(exp.fit$Omega[, 1, 3])),
                        rep('Cor(beta[intercept], beta[size])',
                            2 * length(exp.fit$Omega[, 1, 3])))

tb.cv <- ggplot(baseline.covar, aes(x = value))
tb.cv <- tb.cv + geom_vline(xintercept = 0, colour = 'grey', size = 2)
tb.cv <- tb.cv + geom_histogram(aes(y = ..density..))
tb.cv <- tb.cv + facet_grid(var ~ lab, labeller = label_parsed)
ggsave(tb.cv, filename = '../doc/survival/figure/correlation_marginal.png',
       width = 10, height = 5)

# change in baseline through time
# weibull
wei.grandmean <- wei.fit$mu_prior[, 1]
ww <- list()
for(ii in seq(length(unique(coh)))) {
  ww[[ii]] <- wei.fit$beta[, ii, 1]
}
wei.med <- laply(ww, median)
wei.10 <- laply(ww, function(x) quantile(x, probs = 0.1))
wei.90 <- laply(ww, function(x) quantile(x, probs = 0.9))
wei.bases <- data.frame(stage = seq(length(wei.med)), 
                        med = wei.med, q1 = wei.10, q9 = wei.90)
wei.grand <- data.frame(grand.med = median(wei.grandmean), 
                        grand.q1 = quantile(wei.grandmean, probs = 0.1),
                        grand.q9 = quantile(wei.grandmean, probs = 0.9))
wei.grand$stage <- 1
wei.grand <- rbind(wei.grand, wei.grand)
wei.grand$stage <- c(1, length(wei.med))

# exponential
exp.grandmean <- exp.fit$mu_prior[, 1]
ee <- list()
for(ii in seq(length(unique(coh)))) {
  ee[[ii]] <- exp.fit$beta[, ii, 1]
}
exp.med <- laply(ee, median)
exp.10 <- laply(ee, function(x) quantile(x, probs = 0.1))
exp.90 <- laply(ee, function(x) quantile(x, probs = 0.9))
exp.bases <- data.frame(stage = seq(length(exp.med)), 
                        med = exp.med, q1 = exp.10, q9 = exp.90)
exp.grand <- data.frame(grand.med = median(exp.grandmean), 
                        grand.q1 = quantile(exp.grandmean, probs = 0.1),
                        grand.q9 = quantile(exp.grandmean, probs = 0.9))
exp.grand$stage <- 1
exp.grand <- rbind(exp.grand, exp.grand)
exp.grand$stage <- c(1, length(exp.med))

# put 'em together
wei.bases$label <- 'Weibull'
exp.bases$label <- 'Exponential'
bases <- rbind(wei.bases, exp.bases)
wei.grand$label <- 'Weibull'
exp.grand$label <- 'Exponential'
grand <- rbind(wei.grand, exp.grand)

base.line <- ggplot(bases, aes(x = stage, y = med))
base.line <- base.line + geom_ribbon(data = grand,
                                     mapping = aes(y = grand.med, 
                                                   ymax = grand.q9, 
                                                   ymin = grand.q1), 
                                     alpha = 0.3)
base.line <- base.line + geom_segment(data = grand,
                                      mapping = aes(x = stage[1], 
                                                    xend = stage[2],
                                                    y = grand.med, 
                                                    yend = grand.med))
base.line <- base.line + geom_pointrange(aes(ymax = q9, ymin = q1))
base.line <- base.line + facet_grid(label ~ ., labeller = label_parsed)
base.line <- base.line + labs(x = 'Stage', y = expression(beta[intercept]))
ggsave(base.line, filename = '../doc/survival/figure/intercept_cohort.png',
       width = 10, height = 5)
