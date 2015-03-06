library(rstan)
library(boot)
library(parallel)

source('../R/mung.r')

RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420

# sub sample to test model
library(caret)
dat.full <- dat.full[createDataPartition(dat.full$order, p = 0.1)[[1]], ]

num.occ <- nrow(dat.full)
num.gen <- length(unique(dat.full$genus))
num.ord <- length(unique(as.character(dat.full$order)))

con.gen <- as.numeric(as.factor(dat.full$genus))

gen.ord <- unique(dat.full[, 2:3])
gen.ord <- as.numeric(as.factor(gen.ord[, 2]))
gen.ord <- mapvalues(gen.ord, from = unique(gen.ord), 
                     to = rank(unique(gen.ord)))


data <- list(C = num.occ, G = num.gen, O = num.ord,
             count = dat.full[, 1], genus = con.gen, order = gen.ord,
             off = dat.full$offset)

marine.sample <- stan(file = '../stan/neg_bin_mod.stan')
stan.list <- mclapply(1:4, mc.cores = detectCores(),
                      function(x) stan(fit = marine.sample,
                                       seed = seed,
                                       data = data,
                                       chains = 1, chain_id = x,
                                       refresh = -1))
attempt.one <- sflist2stanfit(stan.list)

#with(data, {stan_rdump(list = c('C', 'G', 'O', 'count', 
#                                'genus', 'order', 'off'),
#                       file = '../data/data_dump/count_info.data.R')})
