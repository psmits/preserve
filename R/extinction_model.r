library(rstan)
library(boot)
library(parallel)

source('../R/mung.r')

# age.data

# remove youngest cohort?

# make cohort, regime, and class into numerics

# stan model
# duration ~ Weibull(alpha, exp(-(intercept + class membership + cohort)
# intercept ~ Normal(0, 10)
#
# class membership ~ Normal(0, sigma)
#
# cohort ~ Normal(regime, sigma)
# regime ~ Normal(0, 10)
