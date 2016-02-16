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
library(ellipse)
library(loo)
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

data <- read_rdump('../data/data_dump/fauna_info.data.R')

pat <- 'faun_'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
#ids <- rep(1:(length(outs) / 4), each = 4)
#outs <- split(outs, ids)

wfit <- read_stan_csv(outs)
#wfit <- llply(outs, read_stan_csv)
#log.liks <- llply(wfit, extract_log_lik)

## this will need to be updated with number of models
#loo.est <- llply(log.liks, loo)
#loo.table <- loo::compare(loo.est[[1]], loo.est[[2]], 
#                          loo.est[[3]], loo.est[[4]])
#
## this will need to be updated with number of models
#waic.est <- llply(log.liks, waic)
#waic.table <- loo::compare(waic.est[[1]], waic.est[[2]], 
#                           waic.est[[3]], waic.est[[4]])
#
#best <- str_extract(names(which.max(waic.table[, ncol(waic.table)])), '[0-9]')
#best <- as.numeric(best)

#wei.fit <- rstan::extract(wfit[[best]], permuted = TRUE)
wei.fit <- rstan::extract(wfit, permuted = TRUE)
wr <- wei.fit$y_tilde[sample(nrow(wei.fit$y_tilde), 1000), ]

# this will need to be updated with number of models
#npred <- ifelse(best %in% c(2, 3), 5, 6)  # 
best <- 3
npred <- 5


#
theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 15),
             legend.text = element_text(size = 15),
             legend.title = element_text(size = 15),
             legend.key.size = unit(1, 'cm'),
             strip.text = element_text(size = 18))

# data setup
coh <- c(data$cohort_unc, data$cohort_cen)
gro <- c(data$group_unc, data$group_cen)
rage <- c(data$occupy_unc, data$occupy_cen)
envs <- c(data$env_unc, data$env_cen)
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
wei.surv <- apply(wr, 1, function(x) survfit(Surv(x) ~ 1))
wei.surv <- llply(wei.surv, function(x) {
                  y <- data.frame(time = x$time, surv = x$surv)
                  y <- rbind(c(0, 1), y[-nrow(y), ])
                  y})
wei.surv <- Reduce(rbind, Map(function(x, y) {
                              x$group <- y
                              x}, 
                              x = wei.surv, 
                              y = seq(length(wei.surv))))
wei.surv$label <- 'Weibull'
# naming/legacy code issue
sim.surv <- wei.surv

# fit model
surv.plot <- ggplot(emp.surv, aes(x = time, y = surv))
surv.plot <- surv.plot + geom_line(data = sim.surv, 
                                   aes(x = time, y = surv, group = group),
                                   colour = 'black', alpha = 0.01)
surv.plot <- surv.plot + geom_line(colour = 'blue')
surv.plot <- surv.plot + coord_cartesian(xlim = c(-0.5, max(duration) + 2))
surv.plot <- surv.plot + labs(x = 'Duration (t)', 
                              y = 'Pr(t < T)')
surv.plot <- surv.plot + theme(axis.title = element_text(size = 25))
surv.plot <- surv.plot + scale_y_continuous(trans=log10_trans(),
                                            breaks = c(0.01, 0.1, 0.5, 1))
ggsave(surv.plot, filename = '../doc/figure/survival_curves.pdf',
       width = 6, height = 5, dpi = 600)

# in b&w
surv.plot <- ggplot(emp.surv, aes(x = time, y = surv))
surv.plot <- surv.plot + geom_line(data = sim.surv, 
                                   aes(x = time, y = surv, group = group),
                                   colour = 'black', alpha = 0.01)
surv.plot <- surv.plot + geom_line(colour = 'grey')
surv.plot <- surv.plot + coord_cartesian(xlim = c(-0.5, max(duration) + 2))
surv.plot <- surv.plot + labs(x = 'Duration (t)', 
                              y = 'Pr(t < T)')
surv.plot <- surv.plot + theme(axis.title = element_text(size = 25))
surv.plot <- surv.plot + scale_y_continuous(trans=log10_trans(),
                                            breaks = c(0.01, 0.1, 0.5, 1))
ggsave(surv.plot, filename = '../doc/figure/survival_curves_bw.pdf',
       width = 6, height = 5, dpi = 600)


# posterior predictive point checks
quant <- laply(wr, function(x) quantile(x, seq(0.1, 0.9, by = 0.05)))
qudur <- quantile(duration, seq(0.1, 0.9, by = 0.05))
qp <- colSums(t(apply(quant, 1, function(x) x > qudur))) / nrow(quant)

point.check <- data.frame(quan = seq(0.1, 0.9, by = 0.05), percent = qp)
point.plot <- ggplot(point.check, aes(x = quan, y = percent))
point.plot <- point.plot + geom_point()
point.plot <- point.plot + labs(x = 'quantile of observed duration',
                                y = 'P(quantile of simulated > 
                                quantile of observed)')
ggsave(point.plot, filename = '../doc/figure/quantile.pdf',
       width = 6, height = 5, dpi = 600)


est.shotgun <- data.frame(obs = duration, 
                          sim = colMeans(wei.fit$hold))

shot.plot <- ggplot(est.shotgun, aes(x = obs, y = sim))
shot.plot <- shot.plot + stat_function(fun = function(x) x, 
                                       lty = 'dashed', 
                                       colour = 'darkgrey')
shot.plot <- shot.plot + geom_point(alpha = 0.5)
shot.plot <- shot.plot + labs(x = 'Observed duration',
                              y = expression(tilde(sigma)))
ggsave(shot.plot, filename = '../doc/figure/shotgun.pdf',
       width = 6, height = 5, dpi = 600)

# quality of fit is medium, though a lot is captured


# make plot of correlation and covariance matrices
# row is sample
# dim 2 is row
# dim 3 is col
get.covcor <- function(stanfit, npred) {
  cor.median <- matrix(, ncol = npred, nrow = npred)
  cor.mean <- matrix(, ncol = npred, nrow = npred)
  cor.10 <- matrix(, ncol = npred, nrow = npred)
  cor.90 <- matrix(, ncol = npred, nrow = npred)
  cov.median <- matrix(, ncol = npred, nrow = npred)
  cov.mean <- matrix(, ncol = npred, nrow = npred)
  cov.10 <- matrix(, ncol = npred, nrow = npred)
  cov.90 <- matrix(, ncol = npred, nrow = npred)
  for(ii in seq(npred)) {
    for(jj in seq(npred)) {
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
               nn <- c('i', 'r', 'e', 'e2', 'm', 's')
               nn <- nn[seq(npred)]
               rownames(x) <- colnames(x) <- nn
               x})
  out
}
wei.covcor <- get.covcor(wei.fit, npred)

# just for the weibull
col1 <- colorRampPalette(c("red", "white", "blue"))
col1<- col1(200)
col2 <- colorRampPalette(c("grey", "white", "grey"))
col2 <- col2(200)

png(file = '../doc/figure/wei_cor_heatmap.png', 
    width = 1500, height = 1500)
my.plotcorr(wei.covcor[[1]], upper.panel = 'number', 
            col = col1[((wei.covcor[[1]] + 1)/2) * 200], 
            mar = rep(0, 4), cex = 4, cex.lab = 4.5,
            npred = npred)
dev.off()

png(file = '../doc/figure/wei_cor_bw.png', 
    width = 1500, height = 1500)
my.plotcorr(wei.covcor[[1]], upper.panel = 'number', 
            col = col2[((wei.covcor[[1]] + 1)/2) * 200], 
            mar = rep(0, 4), cex = 4, cex.lab = 4.5,
            npred = npred)
dev.off()

# effect of all the covariates
# mean cohort effects
efmu <- colMeans(wei.fit$mu_prior)
efmurange <- apply(wei.fit$mu_prior, 2, function(x) 
                   quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9)))

# effects for each cohort
efbeta <- colMeans(wei.fit$beta)
efbetarange <- efbetaprob <- list()
for(ii in seq(data$O)) {
  efbetarange[[ii]] <- apply(wei.fit$beta[, ii, ], 2, function(x) 
                             quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9)))
  efbetaprob[[ii]] <- apply(wei.fit$beta[, ii, ], 2, function(x) 
                            sum(x > 0) / length(x))
}

ef.df <- t(rbind(mean = efmu, efmurange, pred = seq(npred)))
ef.df <- cbind(rbind(ef.df, ef.df), 
               time = c(rep(1, times = npred), rep(data$O, times = npred)))
ef.df <- data.frame(ef.df)
#ef.df$time

efbeta.h <- Map(function(x) t(rbind(mean = efbeta[x, ], efbetarange[[x]])), 
                seq(data$O))
efbeta.h <- Map(function(x) data.frame(efbeta.h[[x]], 
                                       time = x, 
                                       pred = seq(npred)), seq(data$O))
efbeta.df <- Reduce(rbind, efbeta.h)
too <- c('beta[0]', 'beta[r]', 'beta[v]', 'beta[v^2]', 'beta[m]', 
         'beta[s]')[seq(npred)]
efbeta.df$pred <- mapvalues(efbeta.df$pred, 
                            from = seq(npred), 
                            to = too)
ef.df$pred <- mapvalues(ef.df$pred, 
                        from = seq(npred), 
                        to = too)

efbeta.df$pred <- factor(efbeta.df$pred, levels = too)
ef.df$pred <- factor(ef.df$pred, levels = too)

efbeta.df$time <- mapvalues(efbeta.df$time, seq(33), lump[5:(5+33-1), 3])
ef.df$time <- mapvalues(ef.df$time, seq(33), lump[5:(5+33-1), 3])


efbeta.plot <- ggplot(efbeta.df, aes(x = time, y = X50.))
efbeta.plot <- efbeta.plot + geom_pointrange(mapping = aes(ymin = X10., 
                                                           ymax = X90.),
                                             fatten = 2)
efbeta.plot <- efbeta.plot + facet_grid(pred ~ .,
                                        scales = 'free_y', switch = 'y',
                                        labeller = label_parsed)
efbeta.plot <- efbeta.plot + geom_ribbon(data = ef.df, 
                                         mapping = aes(ymin = X10.,
                                                       ymax = X90.),
                                         alpha = 0.2)
efbeta.plot <- efbeta.plot + geom_line(data = ef.df, 
                                       mapping = aes(y = X50.),
                                       alpha = 0.5)
efbeta.plot <- efbeta.plot + labs(x = 'Time (My)', y = 'Parameter estimates')
efbeta.plot <- efbeta.plot + scale_x_reverse()
ggsave(efbeta.plot, filename = '../doc/figure/cohort_series.pdf',
       width = 7.5, height = 10, dpi = 600)



# if best != 1
if(!(best %in% c(1, 2))) {
  #efalrange <- c(mean = mean(wei.fit$alpha_mu), 
  #               quantile(wei.fit$alpha_mu, c(0.1, 0.25, 0.5, 0.75, 0.9)))
  #scalrange <- rbind(mean = colMeans(wei.fit$sigma), 
  #                   apply(wei.fit$sigma, 2, function(x) 
  #                         quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9))))
  #rbind(t(efalrange), t(scalrange))
  #### TODO start here to finish a table with above!

  efalcoh <- rbind(mean = colMeans(wei.fit$alpha_cohort),
                   apply(wei.fit$alpha_cohort, 2, function(x) 
                         quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9))))
  efalcoh <- exp(median(wei.fit$alpha_mu) + efalcoh)

  efalcoh.df <- data.frame(cbind(t(efalcoh), time = seq(data$O)))
  efalcoh.df$time <- mapvalues(efalcoh.df$time, seq(33), lump[5:(5+33-1), 3])

  al.h <- c(exp(quantile(wei.fit$alpha_mu, c(0.1, 0.5, 0.9))))
  al.df <- data.frame(rbind(al.h, al.h), 
                      time = c(max(efalcoh.df$time), min(efalcoh.df$time)))

  efalcoh.plot <- ggplot(efalcoh.df, aes(x = time, y = X50.))
  efalcoh.plot <- efalcoh.plot + geom_pointrange(mapping = aes(ymin = X10., 
                                                               ymax = X90.),
                                                 fatten = 2)
  efalcoh.plot <- efalcoh.plot + geom_ribbon(data = al.df, 
                                             mapping = aes(ymin = X10.,
                                                           ymax = X90.),
                                             alpha = 0.2)
  efalcoh.plot <- efalcoh.plot + geom_line(data = al.df, 
                                           mapping = aes(ymin = X50.),
                                           alpha = 0.5)
  efalcoh.plot <- efalcoh.plot + geom_hline(yintercept = 1, colour = 'darkgrey')
  efalcoh.plot <- efalcoh.plot + scale_x_reverse()
  ggsave(efalcoh.plot, filename = '../doc/figure/shape_series.pdf',
         width = 12.5, height = 10, dpi = 600)
}



# environmental effect
sam <- sample(nrow(wei.fit$lp__), 1000)
quad <- function(x, sam) {
  bet <- wei.fit$mu_prior[sam, 1]
  bet <- bet + (wei.fit$mu_prior[sam, 3] * x) + (wei.fit$mu_prior[sam, 4] * x^2)
  if(!(best %in% 1:2)) {
    -(bet) / exp(wei.fit$alpha_mu[sam])  # depends on if alpha varies by cohort
  } else {
    -(bet) / exp(wei.fit$alpha_trans[sam])  # depends on if alpha varies by cohort
  }
}
quad.mean <- function(x, mcoef) {
  if(!(best %in% 1:2)) {
    -(mcoef[1] + (mcoef[2] * x) + (mcoef[3] * x^2)) / exp(mean(wei.fit$alpha_mu[sam]))
  } else {
    -(mcoef[1] + (mcoef[2] * x) + (mcoef[3] * x^2)) / 
      exp(mean(wei.fit$alpha_trans[sam]))
  }
  # depends on if alpha varies by cohort
}

env.d <- c(data$env_unc, data$env_cen)
val <- seq(from = min(env.d), to = max(env.d), by = 0.01)
quadval <- list()
for(ii in seq(length(sam))) {
  quadval[[ii]] <- data.frame(env = val, resp = quad(val, sam[ii]), sim = ii)
}
quadframe <- Reduce(rbind, quadval)

mcoef <- colMeans(wei.fit$mu_prior[sam, ])[c(1, 3, 4)]
meanquad <- data.frame(env = val, resp = quad.mean(val, mcoef))

# add rug showing observed
#   this addition would overpower the big, by cohort graph
env.obs <- data.frame(env = c(data$env_unc, data$env_cen))

mustache <- ggplot(quadframe, aes(x = env, y = resp, group = sim))
mustache <- mustache + geom_line(alpha = 1 / 100)
mustache <- mustache + geom_line(data = meanquad,
                                 mapping = aes(group = NULL),
                                 colour = 'blue')
mustache <- mustache + geom_rug(data = env.obs,
                                mapping = aes(x = env, y = NULL, group = NULL),
                                sides = 'b', alpha = 0.1)
mustache <- mustache + labs(x = 'Environmental preference (v)', 
                            y = expression(paste('log(', sigma, ')')))
ggsave(mustache, filename = '../doc/figure/env_effect.pdf',
       width = 6, height = 5, dpi = 600)


# by cohort
sam <- sample(nrow(wei.fit$lp__), 100)
bet.coh <- wei.fit$beta[sam, , c(1, 3, 4)]
if(!(best %in% 1:2)) {
  alp.coh <- apply(wei.fit$alpha_cohort[sam, ], 2, function(x) 
                   x + wei.fit$alpha_mu[sam])
} else {
  alp.coh <- wei.fit$alpha_trans[sam]
}
val <- seq(from = min(env.d), to = max(env.d), by = 0.01)
dat <- cbind(1, val, val^2)

coh.est <- list()
for(ii in seq(data$O)) {
  h <- list()
  for(jj in seq(length(val))) {
    # all posterior estimates for env value of dat[1, ]
    if(!(best %in% 1:2)) {
      h[[jj]] <- -(bet.coh[, ii, ] %*% dat[jj, ]) / exp(alp.coh[, ii])
    } else {
      h[[jj]] <- -(bet.coh[, ii, ] %*% dat[jj, ]) / exp(alp.coh[ii])
    }
  }
  coh.est[[ii]] <- h
}

# massage into shape
#   val, resp (V2), sim, coh
stg.name <- as.character(lump[5:(5+33-1), 4])
coh.map <- list()
for(jj in seq(data$O)) {
  h <- Map(function(x, y) {
           cbind(val = x, resp = coh.est[[jj]][[y]], sim = seq(1000))}, 
           x = val, y = seq(length(val)))
  h <- Reduce(rbind, h)
  coh.map[[jj]] <- data.frame(h, coh = stg.name[jj])
}
coh.df <- Reduce(rbind, coh.map)


cohmust <- ggplot(coh.df, aes(x = val, y = V2, group = sim))
cohmust <- cohmust + geom_line(data = meanquad,
                               mapping = aes(x = env,
                                             y = resp,
                                             group = NULL),
                               colour = 'darkgrey')
cohmust <- cohmust + geom_line(alpha = 1 / 100, colour = 'blue')
cohmust <- cohmust + facet_wrap(~ coh, switch = 'x', ncol = 7)
cohmust <- cohmust + theme(axis.text = element_text(size = 6),
                           strip.text = element_text(size = 6))
cohmust <- cohmust + labs(x = 'Environmental preference (v)',
                          y = expression(paste('log(', sigma, ')')))
ggsave(cohmust, filename = '../doc/figure/env_cohort.pdf',
       width = 7.5, height = 8, dpi = 600)
