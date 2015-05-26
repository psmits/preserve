library(rstan)
library(arm)
library(parallel)
source('../R/gts.r')
source('../R/mung.r')

data.file <- list.files('../data', pattern = 'Occs')
fossil <- read.csv(paste0('../data/', data.file))
bibr <- fossil

sight <- space.time(bibr, gts = gts)

num <- nrow(sight)
tim <- ncol(sight)

data <- list(N = num, T = tim, sight = sight)

with(data, {stan_rdump(list = c('N', 'T', 'sight'),
                       file = '../data/data_dump/sight_info.data.R')})
