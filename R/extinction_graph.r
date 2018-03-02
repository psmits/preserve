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

pat <- 'surv_'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)

ids <- rep(1:(length(outs) / 4), each = 4)
outs <- split(outs, ids)

# simple comparison of fits
fits <- llply(outs, read_stan_csv)
log.lik <- llply(fits, extract_log_lik)
loo.est <- llply(log.lik, loo)
waic.est <- llply(log.lik, waic)

# comparison tables
loo.tab <- loo::compare(loo.est[[1]], loo.est[[2]], 
                        loo.est[[3]], loo.est[[4]],
                        loo.est[[5]], loo.est[[6]],
                        loo.est[[7]], loo.est[[8]])
waic.tab <- loo::compare(waic.est[[1]], waic.est[[2]], 
                         waic.est[[3]], waic.est[[4]],
                         waic.est[[5]], waic.est[[6]],
                         waic.est[[7]], waic.est[[8]])
#mn <- c('continuous Weibull w/ interaction', 
#        'continuous Weibull w/o interaction', 
#        'discrete Weibull w/ interaction',
#        'discrete Weibull w/o interaction')
#rownames(loo.tab) <- rownames(waic.tab) <- mn


# move on to the plots
# all the plots for all the models
# allows full comparison between model fits and estimates

# plots for when there is interaction
npred <- 6
# continuous weibull
wei.fit <- rstan::extract(fits[[4]], permuted = TRUE)
posterior.plots(data = data, wei.fit = wei.fit, 
                npred = npred, name = 'cweib_inter')

# discrete weibull
wei.fit <- rstan::extract(fits[[8]], permuted = TRUE)
posterior.plots(data = data, wei.fit = wei.fit, 
                npred = npred, name = 'dweib_inter')

# continuous exponential
wei.fit <- rstan::extract(fits[[2]], permuted = TRUE)
posterior.plots(data = data, wei.fit = wei.fit, 
                npred = npred, name = 'cexp_inter')

# discrete exponential (geometric)
wei.fit <- rstan::extract(fits[[6]], permuted = TRUE)
posterior.plots(data = data, wei.fit = wei.fit, 
                npred = npred, name = 'dexp_inter')





# plots for when there is no interaction
npred <- 5
# continuous weibull
wei.fit <- rstan::extract(fits[[3]], permuted = TRUE)
posterior.plots(data = data, wei.fit = wei.fit, 
                npred = npred, name = 'cweib_base')

# continuous weibull
wei.fit <- rstan::extract(fits[[7]], permuted = TRUE)
posterior.plots(data = data, wei.fit = wei.fit, 
                npred = npred, name = 'dweib_base')

# continuous exponential
wei.fit <- rstan::extract(fits[[1]], permuted = TRUE)
posterior.plots(data = data, wei.fit = wei.fit, 
                npred = npred, name = 'cexp_base')

# discrete exponential (geometric)
wei.fit <- rstan::extract(fits[[5]], permuted = TRUE)
posterior.plots(data = data, wei.fit = wei.fit, 
                npred = npred, name = 'dexp_base')
