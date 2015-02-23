library(rstan)
library(boot)
library(parallel)

source('../R/mung.r')

# age.data
num.samp <- nrow(age.data)
num.orig <- length(unique(age.data$orig))
num.regime <- length(unique(age.data$regime))
num.class <- length(unique(age.data$class))

con.orig <- match(age.data$orig, gts)
con.class <- as.numeric(as.factor(age.data$class))
con.regime <- as.numeric(age.data$regime)


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


stan(file = '../stan/survival_model.stan')

# remove youngest cohort?

# stan model
# duration ~ Weibull(alpha, exp(-(intercept + class membership + cohort)
# intercept ~ Normal(0, 10)
#
# class membership ~ Normal(0, sigma)
#
# cohort ~ Normal(regime, sigma)
# regime ~ Normal(0, 10)
