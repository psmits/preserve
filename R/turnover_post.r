library(plyr)
library(coda)
library(arm)
library(stringr)
library(reshape2)
library(ggplot2)
source('../R/turnover_functions.r')

load('../data/data_dump/occurrence_data.rdata')  # data
load('../data/mcmc_out/turnover_custom.rdata')  # post

# posterior predicitive check
#   if simulate data from same starting point, do we get the same pattern 
#   of **observed** diversity (visually)?
post.check <- replicate(100, posterior.turnover(post = post, data = data), 
                        simplify = FALSE)

true.seen <- obs.seen <- div.true <- div.obs <- list()
for(ii in seq(length(post.check))) {
  true.seen[[ii]] <- lapply(post.check[[ii]], function(x) 
                            check.count(x, set = 'true'))
  obs.seen[[ii]] <- lapply(post.check[[ii]], function(x) 
                           check.count(x, set = 'obs'))

  div.true[[ii]] <- llply(post.check[[ii]], function(x) colSums(x$true))
  div.obs[[ii]] <- llply(post.check[[ii]], function(x) colSums(x$obs))
}
div.dist <- llply(1:4, function(y) 
                  Reduce(rbind, llply(div.obs, function(x) x[[y]])))
div.dist <- Map(function(x) {
                rownames(x) <- seq(nrow(x))
                x}, div.dist)
div.melt <- melt(div.dist)
names(div.melt) <- c('sim', 'year', 'div', 'prov')
div.melt$prov <- factor(div.melt$prov)

obs <- list()
for(jj in seq(data$J)) {
  obs[[jj]] <- cbind(1:data$C, colSums(data$sight[, , jj]), jj)
}
obs <- data.frame(Reduce(rbind, obs))
names(obs) <- c('year', 'div', 'prov')

prov.div <- ggplot(div.melt, aes(x = year, y = div, group = sim))
prov.div <- prov.div + geom_line(alpha = 0.1)
prov.div <- prov.div + geom_line(data = obs, 
                                 mapping = aes(x = year, y = div, group = NULL),
                                 colour = 'blue')
prov.div <- prov.div + facet_grid(prov ~ .)
prov.div <- prov.div + labs(x = 'time', y = 'observed diversity')
ggsave(plot = prov.div, filename = '../doc/gradient/figure/obs_div.pdf',
       width = 10, height = 5)


# "true" diversity...
true.diversity <- function(data, post) {
  samp <- sample(4000, 1)
  regions <- list()
  for(jj in seq(data$J)) {
    hold <- c()
    for(cc in seq(data$C)) {
      hold[cc] <- sum(post$z[[jj]][samp, , cc])
    }
    regions[[jj]] <- hold
  }
  regions
}
est.div <- replicate(1000, true.diversity(data = data, post = post), 
                     simplify = FALSE)
est.div <- llply(seq(data$J), function(y) 
                 Reduce(rbind, llply(est.div, function(x) x[[y]])))

est.div <- Map(function(x) {
               rownames(x) <- seq(nrow(x))
               x}, est.div)
est.div <- melt(est.div)
names(est.div) <- c('sim', 'year', 'div', 'prov')
est.div$prov <- factor(est.div$prov)

prov.est <- ggplot(est.div, aes(x = year, y = div, group = sim))
prov.est <- prov.est + geom_line(alpha = 0.1) + facet_grid(prov ~ .)
prov.est <- prov.est + labs(x = 'time', y = 'estimated diversity')
ggsave(plot = prov.est, filename = '../doc/gradient/figure/est_div.pdf',
       width = 10, height = 5)

# turnover probability
turnover.prob <- function(data, post) {
  samp <- sample(4000, 1)
  regions <- list()
  for(jj in seq(data$J)) {
    hold <- c()
    for(cc in seq(data$C - 1)) {
      hold[cc] <- post$turnover[samp, cc, jj]
    }
    regions[[jj]] <- hold
  }
  regions
}
est.turn <- replicate(1000, turnover.prob(data = data, post = post), 
                     simplify = FALSE)
est.turn <- llply(seq(data$J), function(y) 
                 Reduce(rbind, llply(est.turn, function(x) x[[y]])))

est.turn <- Map(function(x) {
               rownames(x) <- seq(nrow(x))
               x}, est.turn)
est.turn <- melt(est.turn)
names(est.turn) <- c('sim', 'year', 'div', 'prov')
est.turn$prov <- factor(est.turn$prov)

turn.est <- ggplot(est.turn, aes(x = year, y = div, group = sim))
turn.est <- turn.est + geom_line(alpha = 0.1) + facet_grid(prov ~ .)
turn.est <- turn.est + labs(x = 'time', y = 'Pr(z(i,t - 1) = 0 | z(i, t) = 1)')
ggsave(plot = turn.est, filename = '../doc/gradient/figure/turnover.pdf',
       width = 10, height = 5)

# remain/gain probabilities
macro.prob <- function(data, post, ww = 'gamma') {
  samp <- sample(4000, 1)
  regions <- list()
  for(jj in seq(data$J)) {
    hold <- c()
    for(cc in seq(data$C - 1)) {
      hold[cc] <- post[[ww]][samp, cc, jj]
    }
    regions[[jj]] <- hold
  }
  regions
}
est.orig <- replicate(1000, macro.prob(data = data, post = post, ww = 'gamma'), 
                      simplify = FALSE)
est.orig <- llply(seq(data$J), function(y) 
                 Reduce(rbind, llply(est.orig, function(x) x[[y]])))

est.orig <- Map(function(x) {
               rownames(x) <- seq(nrow(x))
               x}, est.orig)
est.orig <- melt(est.orig)
names(est.orig) <- c('sim', 'year', 'div', 'prov')
est.orig$prov <- factor(est.orig$prov)

orig.est <- ggplot(est.orig, aes(x = year, y = div, group = sim))
orig.est <- orig.est + geom_line(alpha = 0.1) + facet_grid(prov ~ .)
orig.est <- orig.est + labs(x = 'time', y = 'Pr(z(i, t) = 1 | z(i, t - 1) = 0)')
ggsave(plot = orig.est, filename = '../doc/gradient/figure/entrance.pdf',
       width = 10, height = 5)

# and the other...
est.surv <- replicate(1000, macro.prob(data = data, post = post, ww = 'gamma'), 
                      simplify = FALSE)

est.surv <- llply(seq(data$J), function(y) 
                 Reduce(rbind, llply(est.surv, function(x) x[[y]])))

est.surv <- Map(function(x) {
               rownames(x) <- seq(nrow(x))
               x}, est.surv)
est.surv <- melt(est.surv)
names(est.surv) <- c('sim', 'year', 'div', 'prov')
est.surv$prov <- factor(est.surv$prov)

surv.est <- ggplot(est.surv, aes(x = year, y = div, group = sim))
surv.est <- surv.est + geom_line(alpha = 0.1) + facet_grid(prov ~ .)
surv.est <- surv.est + labs(x = 'time', y = 'Pr(z(i, t) = 1 | z(i, t - 1) = 1)')
ggsave(plot = surv.est, filename = '../doc/gradient/figure/survival.pdf',
       width = 10, height = 5)



# count the number of gains!
# b = sum (1 - z[i, t - 1]) z[i, t]
birth.count <- function(data, post) {
  samp <- sample(4000, 1)
  birth <- list()
  for(jj in seq(data$J)) {
    hold <- c()
    for(ii in seq(from = 1, to = data$C - 1)) {
      hold[ii] <- sum((1 - post$z[[1]][samp, , ii]) * post$z[[1]][samp, , ii + 1])
    }
    birth[[jj]] <- hold
  }
  birth
}
est.birth <- replicate(1000, birth.count(data = data, post = post), simplify = FALSE)
est.birth <- llply(seq(data$J), function(y) 
                   Reduce(rbind, llply(est.birth, function(x) x[[y]])))

est.birth <- Map(function(x) {
               rownames(x) <- seq(nrow(x))
               x}, est.birth)
est.birth <- melt(est.birth)
names(est.birth) <- c('sim', 'year', 'div', 'prov')
est.birth$prov <- factor(est.birth$prov)

# count the number of losses! 
# d = sum (z[i, t - 1]) (1 - z[i, t])
death.count <- function(data, post) {
  samp <- sample(4000, 1)
  death <- list()
  for(jj in seq(data$J)) {
    hold <- c()
    for(ii in seq(from = 1, to = data$C - 1)) {
      hold[ii] <- sum(post$z[[1]][samp, , ii] * (1 - post$z[[1]][samp, , ii + 1]))
    }
    death[[jj]] <- hold
  }
  death  
}
est.death <- replicate(1000, death.count(data = data, post = post), simplify = FALSE)
est.death <- llply(seq(data$J), function(y) 
                   Reduce(rbind, llply(est.death, function(x) x[[y]])))

est.death <- Map(function(x) {
               rownames(x) <- seq(nrow(x))
               x}, est.death)
est.death <- melt(est.death)
names(est.death) <- c('sim', 'year', 'div', 'prov')
est.death$prov <- factor(est.death$prov)
