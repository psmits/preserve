library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(rstan)
library(survival)
library(stringr)
library(grid)

source('../R/extinction_post_sim.r')

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 20),
             axis.title = element_text(size = 30),
             legend.text = element_text(size = 25),
             legend.title = element_text(size = 26),
             legend.key.size = unit(2, 'cm'),
             strip.text = element_text(size = 20))

# lets make survival curves
duration <- c(data$dur_unc, data$dur_cen)
condition <- c(rep(1, data$N_unc), rep(0, data$N_cen))
condition[duration == 1 & condition == 1] <- 2

emp.surv <- survfit(Surv(time = duration, time2 = duration, 
                         event = condition, type = 'interval') ~ 1)
emp.surv <- data.frame(cbind(time = emp.surv$time, surv = emp.surv$surv))
emp.surv <- rbind(c(0, 1), emp.surv)

sim.surv <- llply(ph, function(x) survfit(Surv(x) ~ 1))
sim.surv <- llply(sim.surv, function(x) {
                  y <- data.frame(time = x$time, surv = x$surv)
                  rbind(c(0, 1), y)
                  y})
sim.surv <- Reduce(rbind, Map(function(x, y) {
                              x$group <- y
                              x}, 
                              x = sim.surv, 
                              y = seq(length(sim.surv))))

surv.plot <- ggplot(emp.surv, aes(x = time, y = surv))
surv.plot <- surv.plot + geom_line(data = sim.surv, 
                                   aes(x = time, y = surv, group = group),
                                   colour = 'grey', alpha = 0.2)
surv.plot <- surv.plot + geom_line(size = 1)
surv.plot <- surv.plot + coord_cartesian(xlim = c(-0.5, max(duration) + 2))
surv.plot <- surv.plot + labs(x = 'Duration in stages', y = 'P(T > t)')

# make plot of correlation and covariance matrices
# row is sample
# dim 2 is row
# dim 3 is col
cor.median <- matrix(, ncol = 4, nrow = 4)
cor.mean <- matrix(, ncol = 4, nrow = 4)
cor.10 <- matrix(, ncol = 4, nrow = 4)
cor.90 <- matrix(, ncol = 4, nrow = 4)
cov.median <- matrix(, ncol = 4, nrow = 4)
cov.mean <- matrix(, ncol = 4, nrow = 4)
cov.10 <- matrix(, ncol = 4, nrow = 4)
cov.90 <- matrix(, ncol = 4, nrow = 4)
for(ii in seq(4)) {
  for(jj in seq(4)) {
    cor.median[ii, jj] <- median(extract.fit$Omega[, ii, jj])
    cor.mean[ii, jj] <- mean(extract.fit$Omega[, ii, jj])
    cor.10[ii, jj] <- quantile(extract.fit$Omega[, ii, jj], probs = .1)
    cor.90[ii, jj] <- quantile(extract.fit$Omega[, ii, jj], probs = .9)
    cov.median[ii, jj] <- median(extract.fit$Sigma[, ii, jj])
    cov.mean[ii, jj] <- mean(extract.fit$Sigma[, ii, jj])
    cov.10[ii, jj] <- quantile(extract.fit$Sigma[, ii, jj], probs = .1)
    cov.90[ii, jj] <- quantile(extract.fit$Sigma[, ii, jj], probs = .9)
  }
}
omega.med <- melt(cor.median)
omega.plot <- ggplot(omega.med, aes(x = Var1, y = Var2, fill = value))
omega.plot <- omega.plot + geom_tile()
omega.plot <- omega.plot + scale_fill_gradient(low = 'white', high = 'black')

sigma.med <- melt(cov.median)
sigma.plot <- ggplot(sigma.med, aes(x = Var1, y = Var2, fill = value))
sigma.plot <- sigma.plot + geom_tile()
sigma.plot <- sigma.plot + scale_fill_gradient2(low = 'blue', 
                                                mid = 'white', 
                                                high = 'red')


trait_baseline.covar <- data.frame(value = extract.fit$Omega[, 1, 3])
tb.cv <- ggplot(trait_baseline.covar, aes(x = value))
tb.cv <- tv.cv + geom_histogram(aes(y = ..density..))

base.grandmean <- extract.fit$mu_prior[, 1]
tt <- list()
for(ii in seq(length(unique(coh)))) {
  tt[[ii]] <- extract.fit$beta[, ii, 1]
}
base.med <- laply(tt, median)
base.10 <- laply(tt, function(x) quantile(x, probs = 0.1))
base.90 <- laply(tt, function(x) quantile(x, probs = 0.9))
bases <- data.frame(stage = seq(length(base.med)), 
                    med = base.med, q1 = base.10, q9 = base.90)
grand <- data.frame(grand.med = median(base.grandmean), 
                    grand.q1 = quantile(base.grandmean, probs = 0.1),
                    grand.q9 = quantile(base.grandmean, probs = 0.9))
grand$stage <- 1
grand <- rbind(grand, grand)
grand$stage <- c(1, length(base.med))
base.line <- ggplot(bases, aes(x = stage, y = med))
base.line <- base.line + geom_ribbon(data = grand,
                                     mapping = aes(y = grand.med, 
                                                   ymax = grand.q9, 
                                                   ymin = grand.q1), 
                                     alpha = 0.3)
base.line <- base.line + geom_segment(data = grand,
                                    mapping = aes(x = stage[1], xend = stage[2],
                                                  y = grand.med, yend = grand.med))
base.line <- base.line + geom_pointrange(aes(ymax = q9, ymin = q1))

#################################
## Exploratory graphs
#
## let's get complicated
## taxon, fauna
#taxon <- c(sepkoski$group_unc, sepkoski$group_cen)
#by.taxon <- split(data.frame(duration, condition), taxon)
#by.taxon <- llply(by.taxon, function(x) survfit(Surv(time = x$duration, 
#                                                     time2 = x$duration,
#                                                     event = x$condition,
#                                                     type = 'interval') ~ 1))
#by.taxon <- llply(by.taxon, function(x) {
#                  oo <- data.frame(time = x$time, surv = x$surv)
#                  oo <- rbind(c(0, 1), oo)
#                  oo})
#by.taxon <- Reduce(rbind, Map(function(x, y) cbind(x, label = rep(y, nrow(x))),
#                              by.taxon, names(by.taxon)))
#
#by.taxon$fauna <- sepkoski$fauna[by.taxon$label]
#
#taxon.surv <- ggplot(by.taxon, aes(x = time, y = surv, colour = label))
#taxon.surv <- taxon.surv + geom_step(direction = 'hv')
#taxon.surv <- taxon.surv + facet_wrap( ~ fauna, ncol = 1)
#
#
## cohort, regime
#temporal <- c(sepkoski$cohort_unc, sepkoski$cohort_cen)
#by.time <- split(data.frame(duration, condition), temporal)
#by.time <- llply(by.time, function(x) survfit(Surv(time = x$duration, 
#                                                   time2 = x$duration,
#                                                   event = x$condition,
#                                                   type = 'interval') ~ 1))
#by.time <- llply(by.time, function(x) {
#                 oo <- data.frame(time = x$time, surv = x$surv)
#                 oo <- rbind(c(0, 1), oo)
#                 oo})
#by.time <- Reduce(rbind, Map(function(x, y) cbind(x, label = rep(y, nrow(x))),
#                             by.time, names(by.time)))
#
#by.time$regime <- sepkoski$regime[by.time$label]
#
#
#time.surv <- ggplot(by.time, aes(x = time, y = surv, colour = label))
#time.surv <- time.surv + geom_step(direction = 'hv')
#time.surv <- time.surv + facet_wrap( ~ regime, ncol = 1)
#
#
#
## explore interactions....
#inter <- interaction(taxon, temporal)
#by.inter <- split(data.frame(duration, condition), inter)
#by.inter <- by.inter[laply(by.inter, function(x) dim(x)[1] > 5)]
#by.inter <- by.inter[laply(by.inter, function(x) any(x[, 2] == 1))]
#by.inter <- llply(by.inter, function(x) survfit(Surv(time = x$duration, 
#                                                   time2 = x$duration,
#                                                   event = x$condition,
#                                                   type = 'interval') ~ 1))
#by.inter <- llply(by.inter, function(x) {
#                 oo <- data.frame(time = x$time, surv = x$surv)
#                 oo <- rbind(c(0, 1), oo)
#                 oo})
#by.inter <- Reduce(rbind, Map(function(x, y) cbind(x, label = rep(y, nrow(x))),
#                             by.inter, names(by.inter)))
#
#by.inter <- cbind(by.inter, 
#                  Reduce(rbind, str_split(as.character(by.inter$label), 
#                                          '\\.')))
#names(by.inter) <- c('time', 'surv', 'label', 'taxon', 'cohort')
#
#by.inter$taxon <- factor(as.numeric(as.character(by.inter$taxon)), 
#                         levels = sort(unique(as.numeric(
#                                  as.character(by.inter$taxon)))))
#by.inter$cohort <- factor(as.numeric(as.character(by.inter$cohort)), 
#                         levels = sort(unique(as.numeric(
#                                  as.character(by.inter$cohort)))))
#by.inter$fauna <- as.factor(sepkoski$fauna[by.inter$taxon])
#by.inter$regime <- as.factor(sepkoski$regime[by.inter$cohort])
#
#time.surv.a <- ggplot(by.inter, aes(x = time, y = surv, 
#                                    group = label))
#time.surv.a <- time.surv.a + geom_step(direction = 'hv', alpha = 0.5)
#time.surv.a <- time.surv.a + facet_grid(regime ~ fauna)
#time.surv.a <- time.surv.a + scale_colour_manual(values = cbp)
#
#time.surv.b <- ggplot(by.inter, aes(x = time, y = surv, 
#                                    group = label))
#time.surv.b <- time.surv.b + geom_step(direction = 'hv', alpha = 0.5)
#time.surv.b <- time.surv.b + facet_grid(fauna ~ regime)
#time.surv.b <- time.surv.b + scale_colour_manual(values = cbp)
