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
source('../R/waic.r')
source('../R/mung.r')
source('../R/multiplot.r')
source('../R/extinction_post_sim.r')
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

pat <- 'faun_weib_[0-9].csv'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
wfit <- read_stan_csv(outs)
wei.fit <- rstan::extract(wfit, permuted = TRUE)

wr <- wei.fit$y_tilde[sample(nrow(wei.fit$y_tilde), 1000), ]

#
theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 15),
             legend.text = element_text(size = 15),
             legend.title = element_text(size = 15),
             legend.key.size = unit(1, 'cm'),
             strip.text = element_text(size = 15))

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
wei.surv <- apply(wr, 1, function(x) survfit(Surv(x) ~ 1))
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
# naming/legacy code issue
sim.surv <- wei.surv

# fit model
surv.plot <- ggplot(emp.surv, aes(x = time, y = surv))
surv.plot <- surv.plot + geom_line(data = sim.surv, 
                                   aes(x = time, y = surv, group = group),
                                   colour = 'black', alpha = 0.01)
surv.plot <- surv.plot + geom_line(size = 0.75, colour = 'blue')
surv.plot <- surv.plot + coord_cartesian(xlim = c(-0.5, max(duration) + 2))
#surv.plot <- surv.plot + facet_grid(. ~ label, labeller = label_parsed)
surv.plot <- surv.plot + labs(x = 'Duration (t)', 
                              y = 'Probability surviving longer than t')
surv.plot <- surv.plot + theme(axis.title = element_text(size = 25),
                               axis.title.y = element_text(size = 20))
ggsave(surv.plot, filename = '../doc/figure/survival_curves.pdf',
       width = 6, height = 5, dpi = 600)

# in b&w
surv.plot <- ggplot(emp.surv, aes(x = time, y = surv))
surv.plot <- surv.plot + geom_line(data = sim.surv, 
                                   aes(x = time, y = surv, group = group),
                                   colour = 'black', alpha = 0.01)
surv.plot <- surv.plot + geom_line(size = 0.75, colour = 'grey')
surv.plot <- surv.plot + coord_cartesian(xlim = c(-0.5, max(duration) + 2))
#surv.plot <- surv.plot + facet_grid(. ~ label, labeller = label_parsed)
surv.plot <- surv.plot + labs(x = 'Duration t', 
                              y = 'Probability surviving greater than t')
surv.plot <- surv.plot + theme(axis.title = element_text(size = 25))
ggsave(surv.plot, filename = '../doc/figure/survival_curves_bw.pdf',
       width = 6, height = 5, dpi = 600)


# posterior predictive point checks
quant <- laply(wr, function(x) quantile(x, seq(0.1, 0.9, by = 0.05)))
qudur <- quantile(duration, seq(0.1, 0.9, by = 0.05))
qp <- colSums(t(apply(quant, 1, function(x) x > qudur))) / nrow(quant)

point.check <- data.frame(quan = seq(0.1, 0.9, by = 0.05), percent = qp)
point.plot <- ggplot(point.check, aes(x = quan, y = percent))
point.plot <- point.plot + geom_point(size = 1.5)
point.plot <- point.plot + labs(x = 'quantile of observed duration',
                                y = 'P(quantile of simulated > 
                                       quantile of observed)')
ggsave(point.plot, filename = '../doc/figure/quantile.pdf',
       width = 6, height = 5, dpi = 600)


est.shotgun <- data.frame(obs = duration, 
                          sim = wei.fit$hold[2, ])

shot.plot <- ggplot(est.shotgun, aes(x = obs, y = sim))
shot.plot <- shot.plot + stat_function(fun = function(x) x, 
                                       size = 1.5,
                                       lty = 'dashed', 
                                       colour = 'grey')
shot.plot <- shot.plot + geom_point(alpha = 0.5, size = 1.5)
shot.plot <- shot.plot + labs(x = 'Observed durations',
                              y = '\\tilde{sigma}')
ggsave(shot.plot, filename = '../doc/figure/shotgun.pdf',
       width = 6, height = 5, dpi = 600)

# quality of fit is medium, though a lot is captured


# make plot of correlation and covariance matrices
# row is sample
# dim 2 is row
# dim 3 is col
get.covcor <- function(stanfit) {
  cor.median <- matrix(, ncol = 8, nrow = 8)
  cor.mean <- matrix(, ncol = 8, nrow = 8)
  cor.10 <- matrix(, ncol = 8, nrow = 8)
  cor.90 <- matrix(, ncol = 8, nrow = 8)
  cov.median <- matrix(, ncol = 8, nrow = 8)
  cov.mean <- matrix(, ncol = 8, nrow = 8)
  cov.10 <- matrix(, ncol = 8, nrow = 8)
  cov.90 <- matrix(, ncol = 8, nrow = 8)
  for(ii in seq(8)) {
    for(jj in seq(8)) {
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
               rownames(x) <- c('i', 'r', 'e', 'e2', 'm', 's', 'sXr', 'sXv')
               colnames(x) <- c('i', 'r', 'e', 'e2', 'm', 's', 'sXr', 'sXv')
               x})
  out
}
wei.covcor <- get.covcor(wei.fit)

# just for the weibull
col1 <- colorRampPalette(c("red", "white", "blue"))
col1<- col1(200)
col2 <- colorRampPalette(c("grey", "white", "grey"))
col2 <- col2(200)

png(file = '../doc/figure/wei_cor_heatmap.png', 
    width = 1500, height = 1500)
my.plotcorr(wei.covcor[[1]], upper.panel = 'number', 
            col = col1[((wei.covcor[[1]] + 1)/2) * 200], 
            mar = rep(0, 4), cex = 4, cex.lab = 4.5)
dev.off()

png(file = '../doc/figure/wei_cor_bw.png', 
    width = 1500, height = 1500)
my.plotcorr(wei.covcor[[1]], upper.panel = 'number', 
            col = col2[((wei.covcor[[1]] + 1)/2) * 200], 
            mar = rep(0, 4), cex = 4, cex.lab = 4.5)
dev.off()

# mean of all coefficients
# sd of all coefficients
param.est <- rbind(data.frame(p = c('mu_i', 'mu_r', 
                                    'mu_v', 'mu_v2', 
                                    'mu_m', 'mu_s',
                                    'mu_rXs', 'mu_vXs'),
                              m = apply(wei.fit$mu_prior, 2, mean), 
                              s = apply(wei.fit$mu_prior, 2, sd),
                              l = apply(wei.fit$mu_prior, 2, 
                                        function(x) quantile(x, 0.1)), 
                              h = apply(wei.fit$mu_prior, 2, 
                                        function(x) quantile(x, 0.9))),
                   data.frame(p = c('tau_i', 'tau_r', 
                                    'tau_v', 'tau_v2', 
                                    'tau_m', 'tau_s',
                                    'tau_rXs', 'tau_vXs'),
                              m = apply(wei.fit$sigma, 2, mean), 
                              s = apply(wei.fit$sigma, 2, sd),
                              l = apply(wei.fit$sigma, 2, 
                                        function(x) quantile(x, 0.1)), 
                              h = apply(wei.fit$sigma, 2, 
                                        function(x) quantile(x, 0.9))))
param.table <- xtable(param.est, label = 'tab:param')
print.xtable(param.table, file = '../doc/table_param.tex')

# histogram of posterior of correlation between inter and env
baseline.covar <- data.frame(value = c(wei.fit$Omega[, 1, 2],
                                       wei.fit$Omega[, 1, 3],
                                       wei.fit$Omega[, 1, 4],
                                       wei.fit$Omega[, 1, 5],
                                       wei.fit$Omega[, 1, 6],
                                       wei.fit$Omega[, 1, 7],
                                       wei.fit$Omega[, 1, 8]),
                             lab = c(rep('Weibull', 
                                         length(wei.fit$Omega[, 1, 3])),
                                     rep('Weibull', 
                                         length(wei.fit$Omega[, 1, 3])),
                                     rep('Weibull', 
                                         length(wei.fit$Omega[, 1, 3])),
                                     rep('Weibull', 
                                         length(wei.fit$Omega[, 1, 3])),
                                     rep('Weibull', 
                                         length(wei.fit$Omega[, 1, 3])),
                                     rep('Weibull', 
                                         length(wei.fit$Omega[, 1, 3])),
                                     rep('Weibull', 
                                         length(wei.fit$Omega[, 1, 3]))))
baseline.covar$var <- c(rep('Cor(beta[0], beta[r])',
                            length(wei.fit$Omega[, 1, 3])),
                        rep('Cor(beta[0], beta[v])',
                            length(wei.fit$Omega[, 1, 3])),
                        rep('Cor(beta[0], beta[v^2])',
                            length(wei.fit$Omega[, 1, 3])),
                        rep('Cor(beta[0], beta[m])',
                            length(wei.fit$Omega[, 1, 3])),
                        rep('Cor(beta[0], beta[s])',
                            length(wei.fit$Omega[, 1, 3])),
                        rep('Cor(beta[0], beta[rXs])',
                            length(wei.fit$Omega[, 1, 3])),
                        rep('Cor(beta[0], beta[vXs])',
                            length(wei.fit$Omega[, 1, 3])))
baseline.covar$var <- factor(baseline.covar$var, 
                             levels = c('Cor(beta[0], beta[r])',
                                        'Cor(beta[0], beta[v])',
                                        'Cor(beta[0], beta[v^2])',
                                        'Cor(beta[0], beta[m])',
                                        'Cor(beta[0], beta[s])',
                                        'Cor(beta[0], beta[rXs])',
                                        'Cor(beta[0], beta[vXs])'))
tb.cv <- ggplot(baseline.covar, aes(x = value))
tb.cv <- tb.cv + geom_vline(xintercept = 0, colour = 'grey', size = 2)
tb.cv <- tb.cv + geom_histogram(aes(y = ..density..))
tb.cv <- tb.cv + facet_grid(var ~ ., labeller = label_parsed)
tb.cv <- tb.cv + labs(x = 'Correlation', y = 'Prob. Density', title = 'B')
tb.cv <- tb.cv + theme(axis.text = element_text(size = 30),
                       axis.title = element_text(size = 40),
                       strip.text = element_text(size = 30),
                       plot.title = element_text(size = 50, hjust = 0))
ggsave(tb.cv, filename = '../doc/figure/correlation_marginal.pdf',
       width = 10, height = 9, dpi = 600)

# mixed figure
png(file = '../doc/figure/cor_mixed.png', 
    width = 3000, height = 1500)
par(mfrow=c(1,2))
my.plotcorr(wei.covcor[[1]], upper.panel = 'number', 
            col = col1[((wei.covcor[[1]] + 1)/2) * 200], 
            #mar = rep(0, 4), 
            cex = 4, cex.lab = 4.5, 
            cex.main = 4.5, main = 'A', adj = 0)
plot.new()
vps <- baseViewports()
pushViewport(vps$figure)
vp1 <- plotViewport(c(1.8, 1, 0, 1))
# plot the ggplot using the print command
print(tb.cv, vp=vp1)
dev.off()

png(file = '../doc/figure/cor_mixed_bw.png', 
    width = 3000, height = 1500)
par(mfrow=c(1,2))
my.plotcorr(wei.covcor[[1]], upper.panel = 'number', 
            col = col2[((wei.covcor[[1]] + 1)/2) * 200], 
            cex = 4, cex.lab = 4.5, 
            cex.main = 4.5, main = 'A', adj = 0)
plot.new()
vps <- baseViewports()
pushViewport(vps$figure)
vp1 <- plotViewport(c(1.8, 1, 0, 1))
# plot the ggplot using the print command
print(tb.cv, vp=vp1)
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

ef.df <- t(rbind(mean = efmu, efmurange, pred = seq(8)))
ef.df <- cbind(rbind(ef.df, ef.df), 
               time = c(rep(1, times = 8), rep(data$O, times = 8)))
ef.df <- data.frame(ef.df)
#ef.df$time

efbeta.h <- Map(function(x) t(rbind(mean = efbeta[x, ], efbetarange[[x]])), 
                seq(data$O))
efbeta.h <- Map(function(x) data.frame(efbeta.h[[x]], 
                                       time = x, 
                                       pred = seq(8)), seq(data$O))
efbeta.df <- Reduce(rbind, efbeta.h)
efbeta.df$pred <- mapvalues(efbeta.df$pred, 
                            from = seq(8), 
                            to = c('beta[0]', 'beta[r]', 
                                   'beta[v]', 'beta[v^2]', 
                                   'beta[m]', 'beta[s]',
                                   'beta[rXs]', 'beta[vXs]'))
ef.df$pred <- mapvalues(ef.df$pred, 
                        from = seq(8), 
                        to = c('beta[0]', 'beta[r]', 
                               'beta[v]', 'beta[v^2]', 
                               'beta[m]', 'beta[s]',
                               'beta[rXs]', 'beta[vXs]'))

efbeta.df$pred <- factor(efbeta.df$pred, 
                         levels = c('beta[0]', 'beta[r]', 
                                    'beta[v]', 'beta[v^2]', 
                                    'beta[m]', 'beta[s]',
                                    'beta[rXs]', 'beta[vXs]'))
ef.df$pred <- factor(ef.df$pred, 
                     levels = c('beta[0]', 'beta[r]', 
                                'beta[v]', 'beta[v^2]', 
                                'beta[m]', 'beta[s]',
                                'beta[rXs]', 'beta[vXs]'))

efbeta.df$time <- mapvalues(efbeta.df$time, seq(33), lump[5:(5+33-1), 3])
ef.df$time <- mapvalues(ef.df$time, seq(33), lump[5:(5+33-1), 3])


efbeta.plot <- ggplot(efbeta.df, aes(x = time, y = X50.))
efbeta.plot <- efbeta.plot + geom_pointrange(mapping = aes(ymin = X10., 
                                                           ymax = X90.),
                                             fatten = 2)
efbeta.plot <- efbeta.plot + facet_grid(pred ~ .,
                                        scales = 'free_y', switch = 'y')
efbeta.plot <- efbeta.plot + geom_ribbon(data = ef.df, 
                                         mapping = aes(ymin = X10.,
                                                       ymax = X90.),
                                         alpha = 0.2)
efbeta.plot <- efbeta.plot + geom_line(data = ef.df, 
                                       mapping = aes(y = X50.),
                                       alpha = 0.5)
efbeta.plot <- efbeta.plot + labs(x = 'Time', y = 'beta')
efbeta.plot <- efbeta.plot + scale_x_reverse()
ggsave(efbeta.plot, filename = '../doc/figure/cohort_series.pdf',
       width = 12.5, height = 15, dpi = 600)




efalrange <- c(mean = mean(wei.fit$alpha_mu), 
               quantile(wei.fit$alpha_mu, c(0.1, 0.25, 0.5, 0.75, 0.9)))
scalrange <- rbind(mean = colMeans(wei.fit$sigma), 
                   apply(wei.fit$sigma, 2, function(x) 
                         quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9))))
rbind(t(efalrange), t(scalrange))
#### TODO start here finish a table!

efalcoh <- rbind(mean = colMeans(wei.fit$alpha_cohort),
                 apply(wei.fit$alpha_cohort, 2, function(x) 
                       quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9))))
efalcoh <- exp(mean(wei.fit$alpha_mu) + efalcoh)

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
efalcoh.plot <- efalcoh.plot + geom_hline(yintercept = 1, colour = 'grey')
efalcoh.plot <- efalcoh.plot + scale_x_reverse()
ggsave(efalcoh.plot, filename = '../doc/figure/shape_series.pdf',
       width = 12.5, height = 10, dpi = 600)



quad <- function(x, sam) {
  bet <- wei.fit$mu_prior[sam, 1]
  bet <- bet + (wei.fit$mu_prior[sam, 3] * x) + (wei.fit$mu_prior[sam, 4] * x^2)
  -(bet) / exp(wei.fit$alpha_mu[sam])
}

val <- seq(from = -0.5, to = 0.5, by = 0.001)
sam <- sample(nrow(wei.fit$alpha), 1000)
quadval <- list()
for(ii in seq(length(sam))) {
  quadval[[ii]] <- data.frame(env = val, resp = quad(val, sam[ii]), sim = ii)
}
quadframe <- Reduce(rbind, quadval)

quad.mean <- function(x, mcoef) {
  -(mcoef[1] + (mcoef[2] * x) + (mcoef[3] * x^2)) / exp(mean(wei.fit$alpha_mu))
}
mcoef <- colMeans(wei.fit$mu_prior)[c(1, 3, 4)]
meanquad <- data.frame(env = val, resp = quad.mean(val, mcoef))

mustache <- ggplot(quadframe, aes(x = env, y = resp, group = sim))
mustache <- mustache + geom_line(alpha = 1 / 100)
mustache <- mustache + geom_line(data = meanquad,
                                 mapping = aes(group = NULL),
                                 colour = 'black', size = 1.5)
mustache <- mustache + labs(x = 'Environmental preference', y = 'log(sigma)')
ggsave(mustache, filename = '../doc/figure/env_effect.pdf',
       width = 6, height = 5, dpi = 600)
