library(plyr)
library(coda)
library(arm)
library(stringr)
source('../R/turnover_functions.r')

load('../data/data_dump/occurrence_data.rdata')  # data
load('../data/mcmc_out/turnover_jags.rdata')  # post.samp

post <- process.coda(post.samp, data)
save(post, file = '../data/mcmc_out/turnover_custom.rdata')
