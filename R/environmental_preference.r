library(rstan)
library(arm)
library(plyr)
library(parallel)
library(ggplot2)
library(reshape2)
library(grid)
source('../R/multiplot.r')

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 5),
             axis.title = element_text(size = 10),
             legend.text = element_text(size = 7.5),
             legend.title = element_text(size = 8),
             legend.key.size = unit(1, 'cm'),
             strip.text = element_text(size = 2.5))

# these are essentially prior predictive tests... 
# what is average occurrence?
# "win" is epicontinental
# "loss" is open ocean

beta.post <- function(n, N, a, b) {
  alp = n + a
  bet = N - n + b + 1
  out <- list(alp, bet)
  out
}

beta.mode <- function(n, N, a, b) {
  alp <- n + a
  bet <- (N - n) + b
  mod <- (alp - 1) / (alp + bet - 2)
  mod
}

set.seed(37)
bck.epi <- 500
bck.off <- 500
bck.theta <- rbeta(bck.epi + bck.off, bck.epi, bck.off)  # background
bck.mode <- beta.mode(bck.epi, bck.epi + bck.off, 1, 1)
samp <- c(1, 2, 3, 4, 5, 10, 25, 50, 100)
nsim <- 1000
prior.a <- 1
prior.b <- 1

sims <- llply(samp, function(x) {
              replicate(n = nsim, expr = {
                        rbinom(1, size = x, sample(bck.theta, 1))})})


# issue with point estimates
# need a lot of occurrences
sim.mode <- Map(function(x, y) beta.mode(x, y, prior.a, prior.b), 
                x = sims, y = samp)
mode.diff <- llply(sim.mode, function(x) x - bck.mode)

# now with full posterior
# when we don't have strong evidence, it is better handeled
sim.post <- Map(function(x, y) {
                pp <- beta.post(x, y, prior.a, prior.b)
                pp <- cbind(pp[[1]], pp[[2]])
                pp},
                x = sims, y = samp)

out <- list()
for(ii in seq(length(samp))) {
  est <- list()
  for(jj in seq(nsim)) {
    post.a <- sim.post[[ii]][jj, 1]
    post.b <- sim.post[[ii]][jj, 2]
    posterior <- qbeta(seq(0, 1, by = 0.001), post.a, post.b)
    bck.post <- qbeta(seq(0, 1, by = 0.001), bck.epi, bck.off)
    est[[jj]] <- posterior - bck.post
  }
  out[[ii]] <- est
}



# compare distributions theta.taxon versus theta.background
# mark modal estimate of theta.taxon
xx <- data.frame(x = seq(0, 1, by = 0.001))
gg.list <- list()
for(ii in seq(length(samp))) {
  o <- sample(nsim, 1)
  post.a <- sim.post[[ii]][o, 1]
  post.b <- sim.post[[ii]][o, 2]
  mod <- sim.mode[[ii]][o]
  gg <- ggplot(xx, aes(x = x)) 
  gg <- gg + stat_function(fun = dbeta, args = list(shape1 = post.a,
                                                    shape2 = post.b), 
                           size = 1, colour = 'goldenrod')
  gg <- gg + stat_function(fun = dbeta, args = list(shape1 = bck.epi, 
                                                    shape2 = bck.off), 
                           size = 1, colour = 'skyblue')
  gg <- gg + geom_vline(xintercept = mod)
  gg <- gg + labs(x = paste(samp[ii]))
  gg.list[[ii]] <- gg
}
png(filename = '../doc/survival/figure/env_post_inspect.png',
    width = 500, height = 1000)
multiplot(plotlist = gg.list, cols = floor(length(gg.list) / 4))
dev.off()

# distribution of modal estimates
diffs.list <- list()
maxet <- 0
for(ii in seq(length(samp))) {
  est.df <- data.frame(x = mode.diff[[ii]])
  maxet <- max(maxet, max(density(est.df$x)$y))
  diffs <- ggplot(est.df, aes(x = x))
  diffs <- diffs + geom_vline(xintercept = 0, colour = 'grey')
  diffs <- diffs + geom_histogram(aes(y = ..density..))
  diffs <- diffs + labs(x = paste(samp[ii]))
  diffs.list[[ii]] <- diffs
}
diffs.list <- llply(diffs.list, function(x) 
                    x + coord_cartesian(xlim = c(-0.6, 0.6), 
                                        ylim = c(-0.5, maxet + 2)))
png(filename = '../doc/survival/figure/env_mode_dist.png',
    width = 500, height = 1000)
multiplot(plotlist = diffs.list, cols = floor(length(samp) / 4))
dev.off()

# for estimates of theta.taxon - theta.background see behavior with sample size
theta.list <- list()
maxet <- 0
for(ii in seq(length(samp))) {
  o <- sample(nsim, 1)
  vline <- mode.diff[[ii]][o]
  est.df <- data.frame(x = out[[ii]][[o]])
  maxet <- max(maxet, max(density(est.df$x)$y))
  theta <- ggplot(est.df, aes(x = x))
  theta <- theta + geom_vline(xintercept = 0, colour = 'grey')
  theta <- theta + geom_vline(xintercept = vline, colour = 'blue')
  theta <- theta + geom_histogram(aes(y = ..density..))
  theta <- theta + labs(x = paste(samp[ii]))
  theta.list[[ii]] <- theta
}
theta.list <- llply(theta.list, function(x) 
                    x + coord_cartesian(xlim = c(-0.6, 0.6), 
                                        ylim = c(-0.5, maxet + 2)))
png(filename = '../doc/survival/figure/env_diff.png',
    width = 500, height = 1000)
multiplot(plotlist = theta.list, cols = floor(length(samp) / 4))
dev.off()
