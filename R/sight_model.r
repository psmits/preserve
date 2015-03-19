library(rstan)
library(boot)
library(parallel)

source('../R/mung.r')

RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420

dat.full <- dat.full[table(dat.full$order) > 1, ]

# sub sample to test model
#library(caret)
#dat.full <- dat.full[createDataPartition(dat.full$order, p = 0.2)[[1]], ]

dat.full$order <- as.character(dat.full$order)
num.occ <- nrow(dat.full)
num.gen <- length(unique(as.character(dat.full$genus)))
num.ord <- length(unique(as.character(dat.full$order)))

con.gen <- as.numeric(as.factor(dat.full$genus))
con.ord <- as.numeric(as.factor(dat.full$order))
counts <- dat.full$count

off <- dat.full$offset

gen.ord <- unique(dat.full[, 2:3])
gen.ord <- as.numeric(as.factor(gen.ord[, 2]))
gen.ord <- mapvalues(gen.ord, 
                     from = unique(gen.ord), 
                     to = rank(unique(gen.ord)))

data <- list(N = num.occ, G = num.gen, O = num.ord,
             count = counts, genus = con.gen, order = gen.ord, off = off)

#marine.sample <- stan(file = '../stan/sampling_model.stan')
#stan.list <- mclapply(1:4, mc.cores = detectCores(),
#                      function(x) stan(fit = marine.sample,
#                                       seed = seed,
#                                       data = data,
#                                       chains = 1, chain_id = x,
#                                       refresh = -1))
#attempt.one <- sflist2stanfit(stan.list)
#summary(attempt.one)[[1]][1:10, 1]

with(data, {stan_rdump(list = c('C', 'G', 'O', 'count', 
                                'genus', 'order', 'off'),
                       file = '../data/data_dump/count_info.data.R')})
