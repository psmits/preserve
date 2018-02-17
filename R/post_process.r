library(ggplot2)
library(survival)
library(arm)
library(rstan)
library(loo)

library(plyr)
library(stringr)
library(reshape2)
library(dplyr)

library(bayesplot)
library(grid)
library(gridBase)
library(gridExtra)
library(ellipse)
library(scales)
source('../R/mung.r')
source('../R/multiplot.r')
source('../R/borrow_plotcorr.r')
source('../R/make_plot.r')
source('../R/stan_utility.R')
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

standata <- read_rdump('../data/data_dump/clean_info.data.R')



pat <- 'clean'
outs <- list.files('../data/mcmc_out', pattern = pat, full.names = TRUE)
fit <- read_stan_csv(outs)
check_all_diagnostics(fit)
post <- rstan::extract(fit, permuted = TRUE)

ppcm <- ppc_stat(standata$y, post$y_tilde, stat = 'mean')
ppcerr <- ppc_error_scatter(standata$y, post$y_tilde[1:6, ])
ppcecdf <- ppc_ecdf_overlay(standata$y, post$y_tilde)
ppcdens <- ppc_dens_overlay(standata$y, post$y_tilde)
ppchist <- ppc_hist(standata$y, post$y_tilde[1:6, ])
ppcsfreq <- ppc_freqpoly(standata$y, post$y_tilde[1:6, ])

# group
ppcmg <- ppc_stat_grouped(standata$y, post$y_tilde, stat = 'mean', group = standata$cohort)
ppcv <- ppc_violin_grouped(standata$y, post$y_tilde, group = standata$cohort)
