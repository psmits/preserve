library(rstan)
library(arm)
library(parallel)

source('../R/mung.r')
sepkoski.data <- sepkoski.data[sepkoski.data$class == 'Rhynchonellata', ]
incl <- sepkoski.data$orig %in% gts[which(gts == 'Changhsingian'):length(gts)]
sepkoski.data <- sepkoski.data[incl, ]

# sepkoski.data
num.samp <- nrow(sepkoski.data)

con.orig <- match(as.character(sepkoski.data$orig), gts)
con.orig <- mapvalues(con.orig, from = unique(con.orig), 
                      to = rank(unique(con.orig)))
num.orig <- length(unique(con.orig))

con.class <- as.numeric(as.factor(sepkoski.data$class))
num.class <- length(unique(sepkoski.data$class))

con.regime <- as.numeric(sepkoski.data$regime) - 1  # kludge
num.regime <- length(unique(sepkoski.data$regime))

con.fauna <- as.numeric(as.factor(sepkoski.data$fauna))
num.fauna <- length(unique(sepkoski.data$fauna))

# beta distribution of genus envrionmental occurrence
beta.a <- sepkoski.data$epi + sepkoski.data$epi.bck
beta.b <- sepkoski.data$off + sepkoski.data$off.bck
ml <- (beta.a - 1) / (beta.a + beta.b - 2)

data <- list(duration = sepkoski.data$duration, group = con.class,
             cohort = con.orig, regime = con.regime, fauna = con.fauna,
             env = rescale(logit(ml)),
             occupy = rescale(logit(sepkoski.data$occupy)))

dead <- sepkoski.data$censored != 1
unc <- llply(data, function(x) x[dead])
cen <- llply(data, function(x) x[!dead])

data <- list(dur_unc = unc$duration,
             group_unc = unc$group,
             fauna_unc = unc$fauna,
             cohort_unc = unc$cohort,
             regime_unc = unc$regime,
             occupy_unc = unc$occupy,
             env_unc = unc$env,
             N_unc = length(unc$duration),
             dur_cen = cen$duration,
             group_cen = cen$group,
             fauna_cen = cen$fauna,
             cohort_cen = cen$cohort,
             regime_cen = cen$regime,
             occupy_cen = cen$occupy,
             env_cen = cen$env,
             N_cen = length(cen$duration))

data$samp_unc <- seq(data$N_unc)
data$samp_cen <- seq(from = data$N_unc + 1,
                     to = data$N_unc + data$N_cen,
                     by = 1)

data$N <- num.samp
data$O <- num.orig
data$R <- num.regime
data$C <- num.class
data$F <- num.fauna
data$regime <- unique(sepkoski.data[, c('orig', 'regime')])[, 2] - 1  # kludge

sepkoski.data$fauna <- as.numeric(as.factor(sepkoski.data$fauna))
data$fauna <- unique(sepkoski.data[, c('class', 'fauna')])[, 2]

with(data, {stan_rdump(list = c('dur_unc', 'group_unc', 'cohort_unc', 
                                'regime_unc', 'N_unc', 'dur_cen', 'group_cen', 
                                'cohort_cen', 'regime_cen', 'N_cen', 'samp_unc', 
                                'samp_cen', 'N', 'O', 'R', 'C', 'F',
                                'occupy_unc', 'occupy_cen',
                                'env_unc', 'env_cen',
                                'fauna', 'regime'),
                       file = '../data/data_dump/fauna_info.data.R')})
