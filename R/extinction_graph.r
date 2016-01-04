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
theme_update(axis.text = element_text(size = 20),
             axis.title = element_text(size = 30),
             legend.text = element_text(size = 25),
             legend.title = element_text(size = 27),
             legend.key.size = unit(1, 'cm'),
             strip.text = element_text(size = 30))

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


# combine both simulation sets
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
#exp.covcor <- get.covcor(exp.fit)

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
                                    'mu_e', 'mu_e2', 'mu_m'),
                              m = apply(wei.fit$mu_prior, 2, mean), 
                              s = apply(wei.fit$mu_prior, 2, sd),
                              l = apply(wei.fit$mu_prior, 2, 
                                        function(x) quantile(x, 0.1)), 
                              h = apply(wei.fit$mu_prior, 2, 
                                        function(x) quantile(x, 0.9))),
                   data.frame(p = c('tau_i', 'tau_r', 
                                    'tau_e', 'tau_e2', 'tau_m'),
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
# put 'em together
wei.bases$label <- 'Weibull'
bases <- wei.bases
wei.grand$label <- 'Weibull'
grand <- wei.grand

renum <- sort(unique(mapvalues(bases$stage, 
                               from = unique(bases$stage), 
                               unique(match(as.character(sepkoski.data$orig),
                                            gts)))))
rename <- gts[renum]
matchy <- match(rename, as.character(lump[, 2]))
rename <- lump[matchy, 3]
rename <- rename - ((diff(lump$age_at_base[c(max(matchy) + 1, matchy)]) / 2) * -1)
bases$stage <- mapvalues(bases$stage, from = unique(bases$stage), rename)
grand$stage <- range(bases$stage)

inter.b <- bases
inter.g <- grand


# then range size
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
# put 'em together
wei.bases$label <- 'Weibull'
bases <- wei.bases
wei.grand$label <- 'Weibull'
grand <- wei.grand

renum <- sort(unique(mapvalues(bases$stage, 
                               from = unique(bases$stage), 
                               unique(match(as.character(sepkoski.data$orig),
                                            gts)))))
rename <- gts[renum]
matchy <- match(rename, as.character(lump[, 2]))
rename <- lump[matchy, 3]
rename <- rename - ((diff(lump$age_at_base[c(max(matchy) + 1, matchy)]) / 2) * -1)
bases$stage <- mapvalues(bases$stage, from = unique(bases$stage), rename)
grand$stage <- range(bases$stage)

rage.b <- bases
rage.g <- grand


# the alpha through time
wei.grandmean <- exp(wei.fit$alpha_mu)
ww <- list()
for(ii in seq(length(unique(coh)))) {
  ww[[ii]] <- wei.fit$alpha[, ii]
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
# put 'em together
wei.bases$label <- 'Weibull'
bases <- wei.bases
wei.grand$label <- 'Weibull'
grand <- wei.grand

renum <- sort(unique(mapvalues(bases$stage, 
                               from = unique(bases$stage), 
                               unique(match(as.character(sepkoski.data$orig),
                                            gts)))))
rename <- gts[renum]
matchy <- match(rename, as.character(lump[, 2]))
rename <- lump[matchy, 3]
rename <- rename - ((diff(lump$age_at_base[c(max(matchy) + 1, matchy)]) / 2) * -1)
bases$stage <- mapvalues(bases$stage, from = unique(bases$stage), rename)
grand$stage <- range(bases$stage)

alph.b <- bases
alph.g <- grand

inter.b$type <- 'beta[intercept]'
rage.b$type <- 'beta[range]'
alph.b$type <- 'alpha'
bases <- rbind(inter.b, rage.b, alph.b)
bases$type <- factor(bases$type, levels = c('beta[intercept]', 'beta[range]', 'alpha'))

inter.g$type <- 'beta[intercept]'
rage.g$type <- 'beta[range]'
alph.g$type <- 'alpha'
grand <- rbind(inter.g, rage.g, alph.g)
grand$type <- factor(grand$type, levels = c('beta[intercept]', 'beta[range]', 'alpha'))


gline <- ggplot(bases, aes(x = stage, y = med))
gline <- gline + geom_ribbon(data = grand,
                                     mapping = aes(y = grand.med, 
                                                   ymax = grand.q9, 
                                                   ymin = grand.q1), 
                                     alpha = 0.3)
gline <- gline + geom_segment(data = grand,
                                      mapping = aes(x = stage[1], 
                                                    xend = stage[2],
                                                    y = grand.med, 
                                                    yend = grand.med))
gline <- gline + geom_pointrange(aes(ymax = q9, ymin = q1))
gline <- gline + scale_x_reverse()
gline <- gline + facet_grid(type ~ ., scales = 'free_y', 
                            labeller = label_parsed)
gline <- gline + labs(x = 'Mya', y = 'Estimate')
gline <- gline + theme(plot.title = element_text(hjust = 0, size = 10),
                       strip.text.y = element_text(angle = 0))
ggsave(gline, filename = '../doc/figure/cohort_series.pdf',
       width = 12.5, height = 10, dpi = 600)


# quadratics plot
sam <- sample(nrow(wei.fit$alpha), 100)
coefs <- data.frame(inter = wei.fit$mu_prior[sam, 1],
                    first = wei.fit$mu_prior[sam, 3], 
                    second = wei.fit$mu_prior[sam, 4],
                    alpha = wei.fit$alpha[sam])
coefplot <- alply(as.matrix(coefs), 1, function(coef) {
                  stat_function(fun = function(x) {
                                coef[1] + coef[2] * x + coef[3] * x^2},
                                colour = 'black',
                                alpha = 0.2)})
lab <- round(sum(coefs$second < 0))
x <- data.frame(x = seq(-5, 5, .0001))
quad <- ggplot(x, aes(x = x)) + coefplot
quad <- quad + stat_function(fun = function(x) {
                             mean(coefs$inter) +
                             mean(coefs$first) * x + 
                             mean(coefs$second) * x^2},
                             colour = 'black',
                             alpha = 1,
                             size = 1.5)
quad <- quad + labs(x = 'environmental preference', 
                    y = 'tilde{sigma}')
ggsave(quad, filename = '../doc/figure/environ_quad.pdf',
       width = 7, height = 5, dpi = 600)

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
