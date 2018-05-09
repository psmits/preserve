library(tidyverse)
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
library(bayesplot)
source('../R/multiplot.r')
source('../R/borrow_plotcorr.r')
source('../R/plot_foo.r')
source('../R/make_plot.r')
source('../R/stan_utility.R')
set.seed(420)

# plotting "rules"
theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 15),
             legend.text = element_text(size = 15),
             legend.title = element_text(size = 15),
             legend.key.size = unit(1, 'cm'),
             strip.text = element_text(size = 15))

# time translation file
lump.file <- list.files('../data', pattern = 'lump')
lump <- read.csv(paste0('../data/', lump.file))
gts <- rev(as.character(lump[, 2]))

# data used to fit the model
data <- read_rdump('../data/data_dump/impute_info.data.R')

pat <- 'surv_cweib_base'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)

# simple comparison of fits
fits <- read_stan_csv(outs)
loglik <- extract_log_lik(fits)
waic1 <- waic(loglik)
loo1 <- loo(loglik)
check_all_diagnostics(fits, max_depth = 15)

# move on to the plots
# all the plots for all the models
# allows full comparison between model fits and estimates

# plots for when there is no interaction
npred <- 5
# continuous weibull
wei.fit <- rstan::extract(fits, permuted = TRUE)
posterior.plots(data = data, wei.fit = wei.fit, 
                npred = npred, lump = lump,
                name = 'cweib_base')


# right censored data form
pat <- 'surv_cweib_cens'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)

# simple comparison of fits
fits <- read_stan_csv(outs)
loglik <- extract_log_lik(fits)
waic2 <- waic(loglik)
loo2 <- loo(loglik)
check_all_diagnostics(fits, max_depth = 15)

# move on to the plots
# all the plots for all the models
# allows full comparison between model fits and estimates

# plots for when there is no interaction
npred <- 5
# continuous weibull
wei.fit <- rstan::extract(fits, permuted = TRUE)

posterior.plots(data = data, wei.fit = wei.fit, 
                npred = npred, lump = lump,
                name = 'cweib_cens', left = TRUE)

tab_waic <- loo::compare(waic1, waic2)
tab_looic <- loo::compare(loo1, loo2)
