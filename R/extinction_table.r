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
ids <- rep(1:(length(outs) / 4), each = 4)
outs <- split(outs, ids)

wfit <- llply(outs, read_stan_csv)
log.liks <- llply(wfit, extract_log_lik)

# this will need to be updated with number of models
loo.est <- llply(log.liks, loo)
loo.table <- loo::compare(loo.est[[1]], #loo.est[[2]], 
                          loo.est[[3]], loo.est[[4]])

# this will need to be updated with number of models
waic.est <- llply(log.liks, waic)
waic.table <- loo::compare(waic.est[[1]], waic.est[[2]], 
                           waic.est[[3]], waic.est[[4]])

best <- str_extract(names(which.max(waic.table[, ncol(waic.table)])), '[0-9]')
best <- as.numeric(best)

wei.fit <- rstan::extract(wfit[[best]], permuted = TRUE)

# this will need to be updated with number of models
npred <- ifelse(best %in% c(2, 3), 5, 8)


# start asking questions and making tables

# mean of all coefficients
# sd of all coefficients
name.mu <- c('mu_i', 'mu_r', 'mu_v', 'mu_v2', 'mu_m', 
             'mu_s', 'mu_rXs', 'mu_vXs')
name.tau <- c('tau_i', 'tau_r', 'tau_v', 'tau_v2', 'tau_m', 
              'tau_s', 'tau_rXs', 'tau_vXs')
qq <- c(0.05, 0.25, 0.5, 0.75, 0.95)
param.est <- rbind(data.frame(p = name.mu[npred],
                              m = apply(wei.fit$mu_prior, 2, mean), 
                              s = apply(wei.fit$mu_prior, 2, sd),
                              l = apply(wei.fit$mu_prior, 2, 
                                        function(x) quantile(x, qq)), 
                              h = apply(wei.fit$mu_prior, 2, 
                                        function(x) quantile(x, qq))),
                   data.frame(p = name.tau[npred],
                              m = apply(wei.fit$sigma, 2, mean), 
                              s = apply(wei.fit$sigma, 2, sd),
                              l = apply(wei.fit$sigma, 2, 
                                        function(x) quantile(x, qq)), 
                              h = apply(wei.fit$sigma, 2, 
                                        function(x) quantile(x, qq))))
param.table <- xtable(param.est, label = 'tab:param')
print.xtable(param.table, file = '../doc/table_param.tex')
