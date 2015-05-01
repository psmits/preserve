library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(arm)
library(rstan)
library(survival)
library(stringr)
library(grid)
source('../R/waic.r')
source('../R/multiplot.r')
source('../R/extinction_post_sim.r')

data <- read_rdump('../data/data_dump/fauna_info.data.R')

pat <- 'faun_weib_[0-9].csv'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
wfit <- read_stan_csv(outs)
wei.fit <- extract(wfit, permuted = TRUE)
weibull.out <- post.sim(data = data, fit = wfit, map = FALSE, expo = FALSE)
wr <- weibull.out[[1]]
wr.res <- weibull.out[[2]]

pat <- 'faun_expo_[0-9].csv'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
efit <- read_stan_csv(outs)
exp.fit <- extract(wfit, permuted = TRUE)
exponential.out <- post.sim(data = data, fit = efit, map = FALSE, expo = TRUE)
er <- exponential.out[[1]]
er.res <- exponential.out[[2]]

#
waic(wfit)$waic < waic(efit)$waic

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 20),
             legend.text = element_text(size = 15),
             legend.title = element_text(size = 16),
             legend.key.size = unit(1, 'cm'),
             strip.text = element_text(size = 7))

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
                  y <- rbind(c(0, 1), y)
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
                  y <- rbind(c(0, 1), y)
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
ggsave(surv.plot, filename = '../doc/survival/figure/survival_curves.png',
       width = 8, height = 5)


# deviance residuals
# change this to be x = duration, y = residual
std.res <- melt(wr.res)
std.res <- std.res[std.res$L1 %in% 1:12, ]
std.res$index <- rep(duration, 12)
res <- ggplot(std.res, aes(x = index, y = value))
res <- res + geom_hline(aes(yintercept = 0), colour = 'grey', size = 1)
res <- res + geom_hline(aes(yintercept = 2), colour = 'grey', size = 1, 
                        linetype = 'dashed')
res <- res + geom_hline(aes(yintercept = -2), colour = 'grey', size = 1, 
                        linetype = 'dashed')
res <- res + geom_point(alpha = 0.5, size = 1, position = 'jitter')
res <- res + facet_wrap( ~ L1, nrow = 3, ncol = 4)
res <- res + labs(x = 'Duration in stages', y = 'Deviance residual')
ggsave(res, filename = '../doc/survival/figure/residual_plot.png',
       width = 8, height = 5)


## posterior predictive point checks
quant <- laply(wr, function(x) quantile(x, seq(0.1, 0.9, by = 0.05)))
qudur <- quantile(duration, seq(0.1, 0.9, by = 0.05))
qp <- colSums(t(apply(quant, 1, function(x) x > qudur))) / nrow(quant)
bad <- which(qp > 0.975 | qp < 0.025)
# quality of fit is weak, though a lot is captured

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
      cor.median[jj, ii] <- median(stanfit$Omega[, jj, ii])
      cor.mean[jj, ii] <- mean(stanfit$Omega[, jj, ii])
      cor.10[jj, ii] <- quantile(stanfit$Omega[, jj, ii], probs = .1)
      cor.90[jj, ii] <- quantile(stanfit$Omega[, jj, ii], probs = .9)
      cov.median[jj, ii] <- median(stanfit$Sigma[, jj, ii])
      cov.mean[jj, ii] <- mean(stanfit$Sigma[, jj, ii])
      cov.10[jj, ii] <- quantile(stanfit$Sigma[, jj, ii], probs = .1)
      cov.90[jj, ii] <- quantile(stanfit$Sigma[, jj, ii], probs = .9)
    }
  }

  out <- list(cor.median, cor.mean, cor.10, cor.90, 
              cov.median, cov.mean, cov.10, cov.90)
  out <- llply(out, function(x) {
               rownames(x) <- c('i', 'r', 'e', 'e2', 'm')
               colnames(x) <- c('i', 'r', 'e', 'e2', 'm')
               x})
  out
}
wei.covcor <- get.covcor(wei.fit)
exp.covcor <- get.covcor(exp.fit)

relab.x <- scale_x_discrete(labels = c('i' = expression(beta[0]), 
                                       'r' = expression(beta[r]),
                                       'e' = expression(beta[v]), 
                                       'e2' = expression(beta[v^2]), 
                                       'm' = expression(beta[m])))
relab.y <- scale_y_discrete(labels = c('i' = expression(beta[intercept]), 
                                       'r' = expression(beta[range]),
                                       'e' = expression(beta[environment]), 
                                       'e2' = expression(beta[environment^2]), 
                                       'm' = expression(beta[size])))
# correlation matrix
omega.med <- melt(list(Exponential = exp.covcor[[1]], 
                       Weibull = wei.covcor[[1]]))
omega.plot <- ggplot(omega.med, aes(x = Var1, y = Var2, fill = value))
omega.plot <- omega.plot + geom_tile()
omega.plot <- omega.plot + geom_text(aes(label = round(value, 2)))
omega.plot <- omega.plot + facet_grid(. ~ L1, labeller = label_parsed)
omega.plot <- omega.plot + scale_fill_gradient2(name = 'Median\nCorrelation',
                                                low = 'blue', 
                                                mid = 'white', 
                                                high = 'red')
omega.plot <- omega.plot + relab.x + relab.y
omega.plot <- omega.plot + labs(x = '', y = '')
omega.plot <- omega.plot + theme(axis.text = element_text(size = 15))
ggsave(omega.plot, filename = '../doc/survival/figure/correlation_heatmap.png',
       width = 10, height = 5)

# just for the weibull
omega.wei <- melt(list(Weibull = wei.covcor[[1]]))
weicor.plot <- ggplot(omega.wei, aes(x = Var1, y = Var2, fill = value))
weicor.plot <- weicor.plot + geom_tile()
weicor.plot <- weicor.plot + geom_text(aes(label = round(value, 2)))
weicor.plot <- weicor.plot + scale_fill_gradient2(name = 'Median\nCorrelation',
                                                  low = 'blue', 
                                                  mid = 'white', 
                                                  high = 'red')
weicor.plot <- weicor.plot + relab.x + relab.y
weicor.plot <- weicor.plot + labs(x = '', y = '')
weicor.plot <- weicor.plot + theme(axis.text = element_text(size = 15))
ggsave(weicor.plot, filename = '../doc/survival/figure/wei_cor_heatmap.png',
       width = 7, height = 5)

# covariance matrix
sigma.med <- melt(list(Exponential = exp.covcor[[5]], 
                       Weibull = wei.covcor[[5]]))
sigma.plot <- ggplot(sigma.med, aes(x = Var1, y = Var2, fill = value))
sigma.plot <- sigma.plot + geom_tile()
sigma.plot <- sigma.plot + geom_text(aes(label = round(value, 2)))
sigma.plot <- sigma.plot + facet_grid(. ~ L1, labeller = label_parsed)
sigma.plot <- sigma.plot + scale_fill_gradient2(name = 'Median\nCovariance', 
                                                low = 'blue', 
                                                mid = 'white', 
                                                high = 'red')
sigma.plot <- sigma.plot + relab.x + relab.y
sigma.plot <- sigma.plot + labs(x = '', y = '')
sigma.plot <- sigma.plot + theme(axis.text = element_text(size = 15))
ggsave(sigma.plot, filename = '../doc/survival/figure/covariance_heatmap.png',
       width = 10, height = 5)


# mean of all coefficients
mid <- c(apply(exp.fit$mu_prior, 2, median),
         apply(wei.fit$mu_prior, 2, median))
top <- c(apply(exp.fit$mu_prior, 2, function(x) quantile(x, probs = 0.9)), 
         apply(wei.fit$mu_prior, 2, function(x) quantile(x, probs = 0.9)))
bot <- c(apply(exp.fit$mu_prior, 2, function(x) quantile(x, probs = 0.1)), 
         apply(wei.fit$mu_prior, 2, function(x) quantile(x, probs = 0.1)))
mids <- data.frame(mid, top, bot)
mids$var <- rep(c('i', 'r', 'e', 'e2', 'm'), 2)
mids$var <- factor(mids$var, levels = unique(mids$var))
mids$mod <- rep(c('Exponential', 'Weibull'), each = 5)

relab.x <- scale_x_discrete(labels = c('i' = expression(mu[intercept]), 
                                       'r' = expression(mu[range]),
                                       'e' = expression(mu[environment]), 
                                       'e2' = expression(mu[environment^2]), 
                                       'm' = expression(mu[size])))

coef.mean <- ggplot(mids, aes(x = var, y = mid))
coef.mean <- coef.mean + geom_hline(aes(yintercept = 0), 
                                    colour = 'grey', size = 2)
coef.mean <- coef.mean + geom_pointrange(aes(ymin = bot, ymax = top, 
                                             colour = mod), size = 0.75,
                                         position = position_jitter(width = 0.1, 
                                                                    height = 0))
coef.mean <- coef.mean + relab.x + labs(x = '', y = expression(hat(mu)))
coef.mean <- coef.mean + scale_colour_manual(values = cbp, name = '')
coef.mean <- coef.mean + theme(axis.text.y = element_text(size = 10),
                               axis.title.y = element_text(angle = 0),
                               axis.text.x = element_text(size = 15),
                               legend.text = element_text(size = 10),
                               legend.position = 'left')
ggsave(coef.mean, filename = '../doc/survival/figure/coef_means.png',
       width = 8, height = 4)


# variance vector for covariates
mid <- c(apply(exp.fit$sigma, 2, median),
         apply(wei.fit$sigma, 2, median))
top <- c(apply(exp.fit$sigma, 2, function(x) quantile(x, probs = 0.9)), 
         apply(wei.fit$sigma, 2, function(x) quantile(x, probs = 0.9)))
bot <- c(apply(exp.fit$sigma, 2, function(x) quantile(x, probs = 0.1)), 
         apply(wei.fit$sigma, 2, function(x) quantile(x, probs = 0.1)))
vars <- data.frame(mid, top, bot)
vars$var <- rep(c('i', 'r', 'e', 'e2', 'm'), 2)
vars$var <- factor(vars$var, levels = unique(vars$var))
vars$mod <- rep(c('Exponential', 'Weibull'), each = 5)

relab.x <- scale_x_discrete(labels = c('i' = expression(tau[intercept]), 
                                       'r' = expression(tau[range]),
                                       'e' = expression(tau[environment]), 
                                       'e2' = expression(tau[environment^2]), 
                                       'm' = expression(tau[size])))

coef.var <- ggplot(vars, aes(x = var, y = mid))
coef.var <- coef.var + geom_hline(aes(yintercept = 0), 
                                  colour = 'grey', size = 2)
coef.var <- coef.var + geom_pointrange(aes(ymin = bot, ymax = top, 
                                           colour = mod), size = 0.75,
                                       position = position_jitter(width = 0.1, 
                                                                  height = 0))
coef.var <- coef.var + relab.x + labs(x = '', y = expression(hat(tau)))
coef.var <- coef.var + scale_colour_manual(values = cbp, name = '')
coef.var <- coef.var + theme(axis.text.y = element_text(size = 10),
                             axis.title.y = element_text(angle = 0),
                             axis.text.x = element_text(size = 15),
                             legend.text = element_text(size = 10),
                             legend.position = 'left')
ggsave(coef.var, filename = '../doc/survival/figure/coef_var.png',
       width = 8, height = 4)


# histogram of posterior of correlation between inter and env
baseline.covar <- data.frame(value = c(exp.fit$Omega[, 1, 2], 
                                       wei.fit$Omega[, 1, 2],
                                       exp.fit$Omega[, 1, 3], 
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
                                         length(wei.fit$Omega[, 1, 3])),
                                     rep('Exponential', 
                                         length(wei.fit$Omega[, 1, 3])),
                                     rep('Weibull', 
                                         length(wei.fit$Omega[, 1, 3]))))
baseline.covar$var <- c(rep('Cor(beta[0], beta[r])',
                            2 * length(exp.fit$Omega[, 1, 3])),
                        rep('Cor(beta[0], beta[v])',
                            2 * length(exp.fit$Omega[, 1, 3])),
                        rep('Cor(beta[0], beta[v^2])',
                            2 * length(exp.fit$Omega[, 1, 3])),
                        rep('Cor(beta[0], beta[m])',
                            2 * length(exp.fit$Omega[, 1, 3])))
baseline.covar$var <- factor(baseline.covar$var, 
                             levels = c('Cor(beta[0], beta[r])',
                                        'Cor(beta[0], beta[v])',
                                        'Cor(beta[0], beta[v^2])',
                                        'Cor(beta[0], beta[m])'))
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
#base.line <- base.line + geom_line()
base.line <- base.line + facet_grid(label ~ ., labeller = label_parsed)
base.line <- base.line + labs(x = 'Stage', y = expression(beta[intercept]))
ggsave(base.line, filename = '../doc/survival/figure/intercept_cohort.png',
       width = 10, height = 5)


# change in range size effect through time
# weibull
wei.grandmean <- wei.fit$mu_prior[, 2]
ww <- list()
for(ii in seq(length(unique(coh)))) {
  ww[[ii]] <- wei.fit$beta[, ii, 2]
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
exp.grandmean <- exp.fit$mu_prior[, 2]
ee <- list()
for(ii in seq(length(unique(coh)))) {
  ee[[ii]] <- exp.fit$beta[, ii, 2]
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

rage.line <- ggplot(bases, aes(x = stage, y = med))
rage.line <- rage.line + geom_ribbon(data = grand,
                                     mapping = aes(y = grand.med, 
                                                   ymax = grand.q9, 
                                                   ymin = grand.q1), 
                                     alpha = 0.3)
rage.line <- rage.line + geom_segment(data = grand,
                                      mapping = aes(x = stage[1], 
                                                    xend = stage[2],
                                                    y = grand.med, 
                                                    yend = grand.med))
rage.line <- rage.line + geom_pointrange(aes(ymax = q9, ymin = q1))
#rage.line <- rage.line + geom_line()
rage.line <- rage.line + facet_grid(label ~ ., labeller = label_parsed)
rage.line <- rage.line + labs(x = 'Stage', y = expression(beta[range]))
ggsave(rage.line, filename = '../doc/survival/figure/range_cohort.png',
       width = 10, height = 5)


# quadratics plot
sam <- sample(nrow(exp.fit$mu_prior), 1000)
coefs <- data.frame(first = wei.fit$mu_prior[sam, 3], 
                    second = wei.fit$mu_prior[sam, 4],
                    alpha = wei.fit$alpha[sam])
coefplot <- alply(as.matrix(coefs), 1, function(coef) {
                  stat_function(fun = function(x) {
                                exp(-(coef[1] * x + coef[2] * x^2) / 
                                    coef[3])},
                                colour = 'grey',
                                alpha = 0.2)})
lab <- round(sum(coefs[, 2] > 0) / nrow(coefs), 2)
x <- data.frame(x = seq(-1, 1, 0.001))
quad <- ggplot(x, aes(x = x)) + coefplot
quad <- quad + stat_function(fun = function(x) {
                             cc <- colMeans(coefs)
                             exp(-(cc[1] * x + cc[2] * x^2) /
                                 cc[3])}, size = 1)
quad <- quad + geom_text(y = 1.75, x = 0, 
                         label = paste(lab), size = 10)
quad <- quad + coord_cartesian(ylim = c(-0.5, 2))
quad <- quad + labs(x = expression(v[i]), 
                    y = expression(paste(tilde(sigma[i])/tilde(sigma))))
ggsave(quad, filename = '../doc/survival/figure/environ_quad.png',
       width = 7, height = 5)

# now do quadratics plot with facet for each cohort
# cohort one
sam <- sample(nrow(exp.fit$mu_prior), 1000)
coef.list <- list()
plotlist <- list()
for(ii in seq(unique(coh))) {
  coefs <- data.frame(first = wei.fit$beta[sam, ii, 3],
                      second = wei.fit$beta[sam, ii, 4],
                      alpha = wei.fit$alpha[sam])
  lab <- round(sum(coefs[, 2] > 0) / nrow(coefs), 2) # probability downward
  cols <- ifelse(mean(coefs[, 2]) < 0, 'red', 'black') # is mean is upward?
  coef.list[[ii]] <- coefs
  ss <- sample(nrow(exp.fit$mu_prior), 100)
  coefs <- data.frame(first = wei.fit$beta[ss, ii, 3],
                      second = wei.fit$beta[ss, ii, 4],
                      alpha = wei.fit$alpha[ss])
  coefplot <- alply(as.matrix(coefs), 1, function(coef) {
                    stat_function(fun = function(x) {
                                  exp(-(coef[1] * x + coef[2] * x^2) / 
                                      coef[3])},
                                  colour = 'grey',
                                  alpha = 0.75)})
  quadcoh <- ggplot(x, aes(x = x)) + coefplot
  quadcoh <- quadcoh + geom_hline(yintercept = 1)
  quadcoh <- quadcoh + geom_text(y = 1.75, x = 0, 
                                 label = paste(lab), size = 10, colour = cols)
  quadcoh <- quadcoh + coord_cartesian(ylim = c(-0.5, 2))
  quadcoh <- quadcoh + labs(x = paste(ii), 
                            y = expression(paste(tilde(sigma[i])/tilde(sigma))))
  plotlist[[ii]] <- quadcoh
}
png(filename = '../doc/survival/figure/cohort_quads.png', 
    width = 3000, height = 1500)
multiplot(plotlist = plotlist, 
          layout = matrix(c(rev(seq(unique(coh))), NA), 
                          ncol = 8, byrow = TRUE))
dev.off()

# do the derivative of the coefficients; get the inflection points
# percent of inflection points greater than 0
#   towards epicontinental
p.epi.best <- laply(coef.list, function(x) 
                    sum((-x[1]) / (x[2] * 2) > 0) / nrow(x))
# which are, on average, up ward facing parabolas
wh.meanworst <- which(laply(coef.list, function(x) mean(x[, 2])) < 0)
wh.midworst <- which(laply(coef.list, function(x) median(x[, 2])) < 0)
# percent of draws with downward facing parabolas
per.best <- laply(coef.list, function(x) sum(x[, 2] > 0) / nrow(x))

# get probability that inflection point isn't in the observed range
#   evidence just looking at one "arm"
#     up or down doesn't actually matter!
#   evidence of approximate linearity?
#   the thing is curved here because it is exponentiated (definition)
#     get around this because i'm working with the log-d coefs
#   maybe just between -0.5 and 0.5?
#     need to look at preferences to see how much is end member
p.linear <- laply(coef.list, function(x) {
                  sum(x[, 1] / (x[, 2] * 2) > 1 | 
                      x[, 1] / (x[, 2] * 2) < -1)}) / nrow(coef.list[[1]])
