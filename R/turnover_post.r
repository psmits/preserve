library(plyr)
library(coda)
library(arm)
library(stringr)
source('../R/turnover_functions.r')
load('../data/data_dump/occurrence_data.rdata')  # data
load('../data/mcmc_out/turnover_jags.rdata')  # post.samp

post <- Reduce(rbind, post.samp)

vars <- unique(llply(str_split(colnames(post), pattern = '\\['), 
                     function(x) x[1]))

new.post <- list()
for(ii in seq(length(vars))) {
  temp <- post[, str_detect(colnames(post), 
                                      paste0(vars[[ii]], '\\['))]
  
  dim.check <- str_count(colnames(temp), pattern = '\\,')

  col.num <- ncol(temp) / data$J
  if(col.num == 1) { 
    new.post[[ii]] <- temp
  } else if (all(dim.check <= 1)){
    holder <- array(dim = c(nrow(temp), col.num, data$J))
    colseq <- list(seq(col.num), seq(from = col.num + 1, to = col.num * data$J))
    for(jj in seq(data$J)) {
      holder[, , jj] <- temp[, colseq[[jj]]]
    }
    new.post[[ii]] <- holder
  } else if (all(dim.check > 1)) {
    holder <- list()
    regionholder <- array(dim = c(nrow(temp), col.num, data$J))
    colseq <- list(seq(col.num), seq(from = col.num + 1, to = col.num * data$J))
    for(jj in seq(data$J)) {
      regionholder <- temp[, colseq[[jj]]]
      
      timeseq <- matrix(ncol = 2, nrow = data$C)
      timeseq[1, ] <- c(1, data$R)
      for(kk in seq(from = 2, to = data$C))
        timeseq[kk, ] <- c((data$R * (kk - 1)) + 1, data$R * kk)
      
      timeholder <- array(dim = c(nrow(temp), data$R, data$C))
      for(cc in seq(from = 1, to = data$C)) {
        timeholder[, , cc] <- regionholder[, timeseq[cc, 1]:timeseq[cc, 2]]
      }
    }
    new.post[[ii]] <- timeholder
  }
}
names(new.post) <- unlist(vars)
post <- new.post


#post.check <- replicate(100, posterior.turnover(post = post, data = data), 
#                        simplify = FALSE)
#
#true.seen <- obs.seen <- div.true <- div.obs <- list()
#for(ii in seq(length(post.check))) {
#  true.seen[[ii]] <- lapply(post.check[[ii]], function(x) 
#                            check.count(x, set = 'true'))
#  obs.seen[[ii]] <- lapply(post.check[[ii]], function(x) 
#                           check.count(x, set = 'observed'))
#
#  div.true[[ii]] <- llply(post.check[[ii]], function(x) colSums(x$true))
#  div.obs[[ii]] <- llply(post.check[[ii]], function(x) colSums(x$obs))
#}
#div.dist <- llply(1:4, function(y) 
#                  Reduce(rbind, llply(div.obs, function(x) x[[y]])))
#
##plot(div.dist[[1]][1, ], type = 'l')
##for(ii in 2:nrow(div.dist[[1]])) lines(div.dist[[1]][ii, ])
##lines(colSums(data$sight[1:data$prov[1], ]), col = 'blue')
#
#
## simulations with even "effort" across provinces
#sim.data = list(P = 4, 
#                prov = seq(from = 1000, to = 4000, by = 1000), 
#                C = data$C,
#                sight = matrix(rep(0, 4000), ncol = 1))
#post.sim <- replicate(100, posterior.turnover(post = post, data = sim.data),
#                        simplify = FALSE)
#true.seen <- obs.seen <- div.true <- div.obs <- list()
#for(ii in seq(length(post.sim))) {
#  true.seen[[ii]] <- lapply(post.sim[[ii]], function(x) 
#                            check.count(x, set = 'true'))
#  obs.seen[[ii]] <- lapply(post.sim[[ii]], function(x) 
#                           check.count(x, set = 'observed'))
#
#  div.true[[ii]] <- llply(post.sim[[ii]], function(x) colSums(x$true))
#  div.obs[[ii]] <- llply(post.sim[[ii]], function(x) colSums(x$obs))
#}
## i want to see the "true" curves
#div.sim <- llply(1:4, function(y) 
#                 Reduce(rbind, llply(div.true, function(x) x[[y]])))
#
##plot(div.dist[[1]][1, ], type = 'l')
##for(ii in 2:nrow(div.dist[[1]])) lines(div.dist[[1]][ii, ])
