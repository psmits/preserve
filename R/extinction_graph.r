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

wr <- wei.fit$y_tilde[sample(4000, 1000), ]

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


## posterior predictive point checks
#quant <- laply(wr, function(x) quantile(x, seq(0.1, 0.9, by = 0.05)))
#qudur <- quantile(duration, seq(0.1, 0.9, by = 0.05))
#qp <- colSums(t(apply(quant, 1, function(x) x > qudur))) / nrow(quant)
#bad <- which(qp > 0.975 | qp < 0.025)
## quality of fit is weak, though a lot is captured


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
                                    'mu_v', 'mu_v2', 'mu_m'),
                              m = apply(wei.fit$mu_prior, 2, mean), 
                              s = apply(wei.fit$mu_prior, 2, sd),
                              l = apply(wei.fit$mu_prior, 2, 
                                        function(x) quantile(x, 0.1)), 
                              h = apply(wei.fit$mu_prior, 2, 
                                        function(x) quantile(x, 0.9))),
                   data.frame(p = c('tau_i', 'tau_r', 
                                    'tau_v', 'tau_v2', 'tau_m'),
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
                                       wei.fit$Omega[, 1, 5]),
                             lab = c(rep('Weibull', 
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
                            length(wei.fit$Omega[, 1, 3])))
baseline.covar$var <- factor(baseline.covar$var, 
                             levels = c('Cor(beta[0], beta[r])',
                                        'Cor(beta[0], beta[v])',
                                        'Cor(beta[0], beta[v^2])',
                                        'Cor(beta[0], beta[m])'))
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

ef.df <- t(rbind(mean = efmu, efmurange, pred = seq(5)))
ef.df <- cbind(rbind(ef.df, ef.df), 
               time = c(rep(1, times = 5), rep(data$O, times = 5)))
ef.df <- data.frame(ef.df)

efbeta.h <- Map(function(x) t(rbind(mean = efbeta[x, ], efbetarange[[x]])), 
                seq(data$O))
efbeta.h <- Map(function(x) data.frame(efbeta.h[[x]], 
                                       time = x, 
                                       pred = seq(5)), seq(data$O))
efbeta.df <- Reduce(rbind, efbeta.h)


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
ggsave(efbeta.plot, filename = '../doc/figure/cohort_series.pdf',
       width = 12.5, height = 10, dpi = 600)


efalrange <- quantile(wei.fit$alpha_mu, c(0.1, 0.25, 0.5, 0.75, 0.9))
efalcoh <- apply(wei.fit$alpha_cohort, 2, function(x) 
                 quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9)))

efalcoh <- apply(efalcoh, 1, function(x) x + efalrange[3])
efalcoh.df <- data.frame(cbind(efalcoh, time = seq(data$O)))

efalcoh.df

efalcoh.plot <- ggplot(efalcoh.df, aes(x = time, y = X50.))
efalcoh.plot <- efalcoh.plot + geom_pointrange(mapping = aes(ymin = X10., 
                                                             ymax = X90.),
                                               fatten = 2)



efsamp <- quantile(wei.fit$gamma, c(0.1, 0.25, 0.5, 0.75, 0.9))

#plot(c(data$samp_unc, data$samp_cen), colMeans(wei.fit$alpha))


## quadratics plot
#sam <- sample(nrow(wei.fit$alpha), 100)
#coefs <- data.frame(inter = wei.fit$mu_prior[sam, 1],
#                    first = wei.fit$mu_prior[sam, 3], 
#                    second = wei.fit$mu_prior[sam, 4],
#                    alpha = wei.fit$alpha[sam])
#coefplot <- alply(as.matrix(coefs), 1, function(coef) {
#                  stat_function(fun = function(x) {
#                                coef[1] + coef[2] * x + coef[3] * x^2},
#                                colour = 'black',
#                                alpha = 0.2)})
#lab <- round(sum(coefs$second < 0))
#x <- data.frame(x = seq(-5, 5, .0001))
#quad <- ggplot(x, aes(x = x)) + coefplot
#quad <- quad + stat_function(fun = function(x) {
#                             mean(coefs$inter) +
#                             mean(coefs$first) * x + 
#                             mean(coefs$second) * x^2},
#                             colour = 'black',
#                             alpha = 1,
#                             size = 1.5)
#quad <- quad + labs(x = 'environmental preference', 
#                    y = 'tilde{sigma}')
#ggsave(quad, filename = '../doc/figure/environ_quad.pdf',
#       width = 7, height = 5, dpi = 600)
#
## now do quadratics plot with facet for each cohort
## cohort one
#renum <- sort(unique(mapvalues(coh, 
#                               from = unique(coh), 
#                               unique(match(as.character(sepkoski.data$orig),
#                                            gts)))))
#rename <- gts[renum]
#rename <- as.character(lump[match(rename, as.character(lump[, 2])), 4])
#for(ii in seq(length(rename))) {
#  rename[ii] <- paste0((length(rename) +1 ) - ii, '. ', rename[ii])
#}
#
#sam <- sample(nrow(wei.fit$mu_prior), 1000)
#x <- data.frame(x = seq(-1, 1, 0.001))
#coef.list <- list()
#plotlist <- list()
#for(ii in seq(unique(coh))) {
#  coefs <- data.frame(first = wei.fit$beta[sam, ii, 3],
#                      second = wei.fit$beta[sam, ii, 4],
#                      alpha = wei.fit$alpha[sam])
#  lab <- round(sum(coefs[, 2] > 0) / nrow(coefs), 2) # probability downward
#  cols <- ifelse(mean(coefs[, 2]) < 0, 'red', 'black') # is mean is upward?
#  coef.list[[ii]] <- coefs
#
#  ss <- sample(nrow(wei.fit$mu_prior), 100)
#  coefs <- data.frame(first = wei.fit$beta[ss, ii, 3],
#                      second = wei.fit$beta[ss, ii, 4],
#                      alpha = wei.fit$alpha[ss, ii])
#  mm <- apply(coefs, 2, median)
#  dd <- apply(coefs, 2, sd)
#  coefplot <- alply(as.matrix(coefs), 1, function(coef) {
#                    stat_function(fun = function(x) {
#                                  exp(-(coef[1] * x + coef[2] * x^2) / 
#                                      coef[3])},
#                                  colour = 'grey',
#                                  alpha = 0.75)})
#  quadcoh <- ggplot(x, aes(x = x)) + coefplot
#  quadcoh <- quadcoh + stat_function(fun = function(x, f, s, a) {
#                                     exp(-(f * x + s * x^2) / a)},
#                                     colour = 'black',
#                                     size = 1,
#                                     args = list(f = mm[1],
#                                                 s = mm[2],
#                                                 a = mm[3]))
#  quadcoh <- quadcoh + stat_function(fun = function(x, f, s, a) {
#                                     exp(-(f * x + s * x^2) / a)},
#                                     colour = 'black', 
#                                     linetype = 'dashed',
#                                     size = 0.7, 
#                                     args = list(f = mm[1] - dd[1],
#                                                 s = mm[2] - dd[2],
#                                                 a = mm[3] - dd[3]))
#  quadcoh <- quadcoh + stat_function(fun = function(x, f, s, a) {
#                                     exp(-(f * x + s * x^2) / a)},
#                                     colour = 'black', 
#                                     linetype = 'dashed',
#                                     size = 0.7, 
#                                     args = list(f = mm[1] + dd[1],
#                                                 s = mm[2] + dd[2],
#                                                 a = mm[3] + dd[3]))
#  quadcoh <- quadcoh + geom_text(y = 1.75, x = 0, 
#                                 label = paste(lab), size = 10)#, colour = cols)
#  quadcoh <- quadcoh + coord_cartesian(ylim = c(0, 2))
#  quadcoh <- quadcoh + labs(x = paste(rename[ii]), 
#                            #y = expression(paste(tilde(sigma[i])/tilde(sigma))))
#                            y = 'mulitplier')
#  plotlist[[ii]] <- quadcoh
#}
#png(file = '../doc//figure/cohort_quads_short.png', 
#    width = 1000, height = 300)
#do.call('grid.arrange', c(rev(plotlist)[25:27], ncol = 3))
#dev.off()
#png(file = '../doc//figure/cohort_quads.png', 
#    width = 3000, height = 1500)
#do.call('grid.arrange', c(rev(plotlist), ncol = 7))
#dev.off()


## do the derivative of the coefficients; get the inflection points
## percent of inflection points greater than 0
##   towards epicontinental
#p.epi.best <- laply(coef.list, function(x) 
#                    sum((-x[1]) / (x[2] * 2) > 0) / nrow(x))
## which are, on average, up ward facing parabolas
#wh.meanworst <- which(laply(coef.list, function(x) mean(x[, 2])) < 0)
#wh.midworst <- which(laply(coef.list, function(x) median(x[, 2])) < 0)
## percent of draws with downward facing parabolas
#per.best <- laply(coef.list, function(x) sum(x[, 2] > 0) / nrow(x))
#
## get probability that inflection point isn't in the observed range
##   evidence just looking at one "arm"
##     up or down doesn't actually matter!
##   evidence of approximate linearity?
##   the thing is curved here because it is exponentiated (definition)
##     get around this because i'm working with the log-d coefs
##   maybe just between -0.5 and 0.5?
##     need to look at preferences to see how much is end member
#p.linear <- laply(coef.list, function(x) {
#                  sum(x[, 1] / (x[, 2] * 2) > 1 | 
#                      x[, 1] / (x[, 2] * 2) < -1)}) / nrow(coef.list[[1]])
#save(wei.fit,
#     p.epi.best, 
#     wh.meanworst, 
#     wh.midworst, 
#     per.best, 
#     p.linear, 
#     file = '../data/epi_over_off.rdata')
