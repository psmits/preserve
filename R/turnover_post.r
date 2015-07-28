library(plyr)
library(coda)
library(arm)
library(stringr)
library(reshape2)
library(ggplot2)
source('../R/turnover_functions.r')

load('../data/data_dump/occurrence_data.rdata')  # data
load('../data/mcmc_out/turnover_custom.rdata')  # post


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
div.dist <- llply(1:2, function(y) 
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
prov.div <- prov.div + geom_line() + facet_grid(prov ~ .)
prov.div <- prov.div + geom_line(data = obs, 
                                 mapping = aes(x = year, y = div, group = NULL),
                                 colour = 'blue')
