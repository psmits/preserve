library(rstan)
library(boot)
library(parallel)

source('../R/mung.r')

num.occ <- nrow(dat.full)
num.gen <- length(unique(dat.full$genus))
num.ord <- length(unique(dat.full$order))

con.gen <- as.numeric(as.factor(dat.full$genus))

gen.ord <- unique(dat.full[, 2:3])
gen.ord <- as.numeric(as.factor(gen.ord[, 2]))

data <- list(C = num.occ, G = num.gen, O = num.ord,
             count = dat.full[, 1], genus = con.gen, order = gen.ord,
             off = dat.full$offset)

with(data, {stan_rdump(list = c('C', 'G', 'O', 'count', 
                                'genus', 'order', 'off'),
                       file = '../data/data_dump/count_info.data.R')})
