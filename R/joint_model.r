library(rstanarm)


mod1 <- stan_jm(formulaLong = logBili ~ year + (year | id), 
                dataLong = pbcLong,
                formulaEvent = survival::Surv(futimeYears, death) ~ sex + trt, 
                dataEvent = pbcSurv,
                assoc = c('etavalue', 'etaslope'),
                time_var = "year",
                chains = 1, refresh = 2000, seed = 12345)


