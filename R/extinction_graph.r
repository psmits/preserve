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
source('../R/make_plot.r')
set.seed(420)

# plotting "rules"
theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 15),
             axis.title = element_text(size = 20),
             legend.text = element_text(size = 20),
             legend.title = element_text(size = 20),
             legend.key.size = unit(1, 'cm'),
             strip.text = element_text(size = 20))

# bring in data files and get them prepped
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

data <- read_rdump('../data/data_dump/impute_info.data.R')

pat <- 'faun_'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)

ids <- rep(1:(length(outs) / 4), each = 4)
outs <- split(outs, ids)

#log.liks <- llply(wfit, extract_log_lik)
## this will need to be updated with number of models
#loo.est <- llply(log.liks, loo)
#loo.table <- loo::compare(loo.est[[1]], loo.est[[2]], loo.est[[3]])
#
## this will need to be updated with number of models
#waic.est <- llply(log.liks, waic)
#waic.table <- loo::compare(waic.est[[1]], waic.est[[2]], waic.est[[3]])


# move on to the plots

# plots for when there is interaction
npred <- 6
# continuous weibull
wfit <- read_stan_csv(outs[[3]])
wei.fit <- rstan::extract(wfit, permuted = TRUE)
posterior.plots(data = data, wei.fit = wei.fit, npred = npred, name = 'cweib')

# discrete weibull
wfit <- read_stan_csv(outs[[1]])
wei.fit <- rstan::extract(wfit, permuted = TRUE)
posterior.plots(data = data, wei.fit = wei.fit, npred = npred, name = 'dweib')


## plots for when there is no interaction
#npred <- 5
## continuous weibull
#wfit <- read_stan_csv(outs[[3]])
#wei.fit <- rstan::extract(wfit, permuted = TRUE)
#posterior.plots(data = data, wei.fit = wei.fit, npred = npred, name = 'cweib')
#
## discrete weibull
#wfit <- read_stan_csv(outs[[3]])
#wei.fit <- rstan::extract(wfit, permuted = TRUE)
#posterior.plots(data = data, wei.fit = wei.fit, npred = npred, name = 'dweib')
