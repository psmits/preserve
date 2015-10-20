library(plyr)
library(coda)
library(arm)
library(stringr)
library(reshape2)
library(ggplot2)
library(igraph)
library(grid)
library(scales)
source('../R/turnover_functions.r')
source('../R/turnover_plot_foo.r')
set.seed(420)

load('../data/data_dump/occurrence_data.rdata')  # data
load('../data/mcmc_out/turnover_custom.rdata')  # post

lump.file <- list.files('../data', pattern = 'lump')
lump <- read.csv(paste0('../data/', lump.file))
#
theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 20),
             legend.text = element_text(size = 15),
             legend.title = element_text(size = 16),
             legend.key.size = unit(1, 'cm'),
             strip.text = element_text(size = 15))


# posterior predicitive check
#   if simulate data from same starting point, do we get the same pattern 
#   of **observed** diversity (visually)?
post.check <- replicate(1000, posterior.turnover(post = post, data = data), 
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
for(jj in seq(data$nprov)) {
  obs[[jj]] <- cbind(1:data$nyear, colSums(data$y[, , jj]), jj)
}
obs <- data.frame(Reduce(rbind, obs))
names(obs) <- c('year', 'div', 'prov')

# province names
div.melt$prov <- mapvalues(div.melt$prov, unique(div.melt$prov), 
                           c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
div.melt$prov <- factor(div.melt$prov, levels = 
                        c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
obs$prov <- mapvalues(obs$prov, unique(obs$prov), 
                      c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
obs$prov <- factor(obs$prov, levels = 
                   c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))

# make in years
time.slice <- lump[seq(from = 5, to = data$nyear + 4), ]
div.melt$year <- mapvalues(div.melt$year, 
                           unique(div.melt$year), time.slice[, 3])
obs$year <- mapvalues(obs$year, unique(obs$year), time.slice[, 3])

prov.div <- ggplot(div.melt, aes(x = year, y = div, group = sim))
prov.div <- prov.div + geom_line(alpha = 0.01)
prov.div <- prov.div + geom_line(data = obs, 
                                 mapping = aes(x = year, y = div, group = NULL),
                                 colour = 'blue')
prov.div <- prov.div + facet_grid(prov ~ .)
prov.div <- prov.div + scale_x_reverse()
prov.div <- prov.div + coord_cartesian(ylim = c(0, 250))
prov.div <- prov.div + labs(x = 'time', y = 'observed diversity')
ggsave(plot = prov.div, filename = '../doc/gradient/figure/obs_div.png',
       width = 10, height = 5)

# "true" diversity...
div.dist <- llply(1:4, function(y) 
                  Reduce(rbind, llply(div.true, function(x) x[[y]])))
div.dist <- Map(function(x) {
                rownames(x) <- seq(nrow(x))
                x}, div.dist)
div.mean <- melt(laply(div.dist, colMeans))
names(div.mean) <- c('prov', 'year', 'div')

div.melt <- melt(div.dist)
names(div.melt) <- c('sim', 'year', 'div', 'prov')
div.melt$prov <- factor(div.melt$prov)

# province names
div.melt$prov <- mapvalues(div.melt$prov, unique(div.melt$prov), 
                           c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
div.melt$prov <- factor(div.melt$prov, levels = 
                        c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
div.mean$prov <- mapvalues(div.mean$prov, unique(div.mean$prov), 
                           c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
div.mean$prov <- factor(div.mean$prov, levels = 
                        c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))

# make in years
time.slice <- lump[seq(from = 5, to = data$nyear + 4), ]
div.melt$year <- mapvalues(div.melt$year, 
                           unique(div.melt$year), time.slice[, 3])
div.mean$year <- mapvalues(div.mean$year, 
                           unique(div.mean$year), time.slice[, 3])

prov.est <- ggplot(div.melt, aes(x = year, y = div, group = sim))
prov.est <- prov.est + geom_line(alpha = 0.01)
prov.est <- prov.est + geom_line(data = div.mean,
                                 mapping = aes(x = year, y = div, group = NULL), 
                                 colour = 'blue')
prov.est <- prov.est + scale_x_reverse() + facet_grid(prov ~ .)
prov.est <- prov.est + coord_cartesian(ylim = c(0, 250))
prov.est <- prov.est + labs(x = 'time', y = 'estimated diversity')
ggsave(plot = prov.est, filename = '../doc/gradient/figure/true_div.png',
       width = 10, height = 5)


# need to make an easier to understand figure
#   relative diversity -- lets identification of "gradient"
div.dist <- llply(1:4, function(y) 
                  Reduce(rbind, llply(div.true, function(x) x[[y]])))
div.dist <- Map(function(x) {
                rownames(x) <- seq(nrow(x))
                x}, div.dist)
mean.div <- laply(div.dist, colMeans)
mean.sums <- colSums(mean.div)
for(ii in seq(ncol(mean.div))) {
  mean.div[, ii] <- mean.div[, ii] / mean.sums[ii]
}

time.slice <- lump[seq(from = 5, to = data$nyear + 5), ]
width <- diff(time.slice[, 3]) * -1
time.slice <- time.slice[-nrow(time.slice), ]

colnames(mean.div) <- time.slice[, 3]
melt.mean <- melt(mean.div)

melt.mean$width = 0
for(ii in seq(nrow(time.slice))) {
  melt.mean$width[melt.mean[, 2] == time.slice[ii, 3]] <- width[ii] * 2
}
melt.mean$Var1 <- mapvalues(melt.mean$Var1, 
                            unique(melt.mean$Var1), 
                            c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
melt.mean$Var1 <- factor(melt.mean$Var1, 
                         levels = c('S. Temp', 'S. Trop', 'N. Trop', 'N. Temp'))

rel.plot <- ggplot(melt.mean, aes(x = Var2, y = Var1, fill = value, width = width))
rel.plot <- rel.plot + geom_tile()
rel.plot <- rel.plot + scale_fill_gradient(low = 'white', high = 'black')
rel.plot <- rel.plot + scale_x_reverse() + labs(x = 'Time (Mya)', y = 'Region')
ggsave(plot = rel.plot, filename = '../doc/gradient/figure/rel_div.png',
       width = 10, height = 5)


# occupancy probability??
#   uses psi
#
est.occ <- replicate(1000, occ.prob(post = post, data = data), simplify = FALSE)
est.growth <- growth.rate(est.occ, data)
macro.plot(posterior = est.growth, 
           time = lump, 
           label = 'Growth rate', 
           filename = '../doc/gradient/figure/growth.png',
           process = FALSE,
           log = TRUE)

macro.plot(posterior = est.occ, 
           time = lump, 
           label = 'Occupancy probability', 
           filename = '../doc/gradient/figure/occupancy.png')

# macro probabilities
est.turn <- replicate(1000, macro.prob(data = data, 
                                       post = post, 
                                       ww = 'turnover'), 
                     simplify = FALSE)
macro.plot(posterior = est.turn, time = lump)

est.orig <- replicate(1000, macro.prob(data = data, 
                                       post = post, 
                                       ww = 'gamma'), 
                      simplify = FALSE)
macro.plot(posterior = est.orig, 
           time = lump, 
           label = 'Pr(z(i, t) = 1 | z(i, t - 1) = 0)',
           filename = '../doc/gradient/figure/entrance.png')

est.surv <- replicate(1000, macro.prob(data = data, 
                                       post = post, 
                                       ww = 'phi'), 
                      simplify = FALSE)
macro.plot(posterior = est.surv, 
           time = lump, 
           label = 'Pr(z(i, t) = 0 | z(i, t - 1) = 1)',
           filename = '../doc/gradient/figure/survival.png',
           subtract = TRUE)

est.obs <- replicate(1000, macro.prob(data = data, 
                                      post = post, 
                                      ww = 'p'), 
                      simplify = FALSE)
macro.plot(posterior = est.obs, 
           time = lump, 
           label = 'Pr(y(i, t) = 1 | z(i, t) = 1)',
           '../doc/gradient/figure/observation.png')


# div dep graphs
# beta is orig
orig.test <- orig.sd <- orig.eff <- matrix(ncol = data$nprov, nrow = data$nprov)
for(ii in seq(data$nprov)) {
  orig.eff[ii, ] <- colMeans(post$beta[, , ii])
  orig.sd[ii, ] <- apply(post$beta[, , ii], 2, sd)

  for(jj in seq(data$nprov)) {
    if(orig.eff[ii, jj] > 0) {
      orig.test[ii, jj] <- orig.eff[ii, jj] - 1 * orig.sd[ii, jj]
    } else {
      orig.test[ii, jj] <- orig.eff[ii, jj] + 1 * orig.sd[ii, jj]
    }
  }
}
orig.eff[sign(orig.eff) != sign(orig.test)] = 0

orig.graph <- graph_from_adjacency_matrix(orig.eff, 
                                          mode = 'directed', 
                                          weighted = TRUE)
orig.col <- ifelse(E(orig.graph)$weight <= 0, 'red', 'blue')

#png(file = '../doc/gradient/figure/orig_graph.png',
#    width = 300, height = 150)
plot(orig.graph, 
     edge.width = abs(E(orig.graph)$weight) * 100, 
     edge.color = orig.col,
     edge.curved = TRUE, layout = layout_in_circle(orig.graph))
#dev.off()

# alpha is surv
exti.test <- exti.sd <- exti.eff <- matrix(ncol = data$nprov, nrow = data$nprov)
for(ii in seq(data$nprov)) {
  exti.eff[ii, ] <- colMeans(post$alpha[, , ii])
  exti.sd[ii, ] <- apply(post$alpha[, , ii], 2, sd)

  for(jj in seq(data$nprov)) {
    if(exti.eff[ii, jj] > 0) {
      exti.test[ii, jj] <- exti.eff[ii, jj] - 1 * exti.sd[ii, jj]
    } else {
      exti.test[ii, jj] <- exti.eff[ii, jj] + 1 * exti.sd[ii, jj]
    }
  }
}
exti.eff[sign(exti.eff) != sign(exti.test)] = 0

exti.graph <- graph_from_adjacency_matrix(exti.eff, 
                                          mode = 'directed', 
                                          weighted = TRUE)
exti.col <- ifelse(E(exti.graph)$weight <= 0, 'red', 'blue')

#png(file = '../doc/gradient/figure/exti_graph.png',
#    width = 300, height = 150)
plot(exti.graph, 
     edge.width = abs(E(exti.graph)$weight) * 100,
     edge.color = exti.col,
     edge.curved = TRUE, layout = layout_in_circle(exti.graph))
#dev.off()


## make a rate
#width <- rev(diff(rev(lump[seq(from = 5, to = data$nyear + 4), 3])))
#split.surv <- split(est.surv, est.surv$year)
#est.surv$rate <- 0
#for(ii in seq(length(split.surv))) {
#  ss <- split.orig[[ii]]$div
#  est.surv$rate <- -log(1 - ss) / width[ii]
#}
#
#surv.est <- ggplot(est.surv, aes(x = year, y = rate, group = sim))
#surv.est <- surv.est + geom_line(alpha = 0.01) + facet_grid(prov ~ .)
#surv.est <- surv.est + scale_x_reverse() + coord_trans(y = 'log')
#surv.est <- surv.est + labs(x = 'time', y = 'Extinction/loss rate')
#ggsave(plot = surv.est, filename = '../doc/gradient/figure/surv_rate.png',
#       width = 10, height = 5)



## count the number of gains!
## b = sum (1 - z[i, t - 1]) z[i, t]
#birth.count <- function(data, post) {
#  samp <- sample(4000, 1)
#  birth <- list()
#  for(jj in seq(data$nprov)) {
#    hold <- c()
#    for(ii in seq(from = 1, to = data$nyear - 1)) {
#      hold[ii] <- sum((1 - post$z[[1]][samp, , ii]) * 
#                      post$z[[1]][samp, , ii + 1])
#    }
#    birth[[jj]] <- hold
#  }
#  birth
#}
#est.birth <- replicate(1000, birth.count(data = data, post = post), 
#                       simplify = FALSE)
#est.birth <- llply(seq(data$nprov), function(y) 
#                   Reduce(rbind, llply(est.birth, function(x) x[[y]])))
#
#est.birth <- Map(function(x) {
#               rownames(x) <- seq(nrow(x))
#               x}, est.birth)
#est.birth <- melt(est.birth)
#names(est.birth) <- c('sim', 'year', 'div', 'prov')
#est.birth$prov <- factor(est.birth$prov)
#
#
## count the number of losses! 
## d = sum (z[i, t - 1]) (1 - z[i, t])
#death.count <- function(data, post) {
#  samp <- sample(4000, 1)
#  death <- list()
#  for(jj in seq(data$nprov)) {
#    hold <- c()
#    for(ii in seq(from = 1, to = data$nyear - 1)) {
#      hold[ii] <- sum(post$z[[1]][samp, , ii] * 
#                      (1 - post$z[[1]][samp, , ii + 1]))
#    }
#    death[[jj]] <- hold
#  }
#  death  
#}
#est.death <- replicate(1000, death.count(data = data, post = post), 
#                       simplify = FALSE)
#est.death <- llply(seq(data$nprov), function(y) 
#                   Reduce(rbind, llply(est.death, function(x) x[[y]])))
#
#est.death <- Map(function(x) {
#               rownames(x) <- seq(nrow(x))
#               x}, est.death)
#est.death <- melt(est.death)
#names(est.death) <- c('sim', 'year', 'div', 'prov')
#est.death$prov <- factor(est.death$prov)
