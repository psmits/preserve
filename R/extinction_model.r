library(rstan)
library(boot)
library(parallel)

source('../R/mung.r')

# sub sample to test model
#library(caret)
#age.data <- age.data[createDataPartition(age.data$orig, p = 0.25)[[1]], ]

# age.data
num.samp <- nrow(age.data)

con.orig <- match(as.character(age.data$orig), gts)
con.orig <- mapvalues(con.orig, from = unique(con.orig), 
                      to = rank(unique(con.orig)))
num.orig <- length(unique(con.orig))

con.class <- as.numeric(as.factor(age.data$class))
num.class <- length(unique(age.data$class))

con.regime <- as.numeric(age.data$regime)
num.regime <- length(unique(age.data$regime))


data <- list(duration = age.data$duration, group = con.class,
             cohort = con.orig, regime = con.regime)

dead <- age.data$censored != 1
unc <- llply(data, function(x) x[dead])
cen <- llply(data, function(x) x[!dead])

data <- list(dur_unc = unc$duration,
             group_unc = unc$group,
             cohort_unc = unc$cohort,
             regime_unc = unc$regime,
             N_unc = length(unc$duration),
             dur_cen = cen$duration,
             group_cen = cen$group,
             cohort_cen = cen$cohort,
             regime_cen = cen$regime,
             N_cen = length(cen$duration))

data$samp_unc <- seq(data$N_unc)
data$samp_cen <- seq(from = data$N_unc + 1,
                     to = data$N_unc + data$N_cen,
                     by = 1)

data$N <- num.samp
data$O <- num.orig
data$R <- num.regime
data$C <- num.class
data$regime <- unique(age.data[, c('orig', 'regime')])[, 2]

with(data, {stan_rdump(list = c('dur_unc', 'group_unc', 'cohort_unc', 
                                'regime_unc', 'N_unc', 'dur_cen', 'group_cen', 
                                'cohort_cen', 'regime_cen', 'N_cen', 'samp_unc', 
                                'samp_cen', 'N', 'O', 'R', 'C', 'regime'),
                       file = '../data/data_dump/survival_info.data.R')})
#marine.survival <- stan(file = '../stan/survival_model.stan')
#attempt.one <- stan(fit = marine.survival, data = data)
