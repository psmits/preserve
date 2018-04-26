library(tidyverse)
library(rstan)
source('../R/gts.r')
source('../R/jade.r')
source('../R/mung.r')

data.file <- list.files('../data', pattern = 'Occs')
fossil <- read_csv(paste0('../data/', data.file))
bibr <- fossil

lump.file <- list.files('../data', pattern = 'lump')
lump <- read_csv(paste0('../data/', lump.file))

payne <- read_tsv('../data/payne_bodysize/Occurrence_PaleoDB.txt')

taxon = c('Rhynchonellata', 'Strophomenata', 'Chileata', 'Obolellata', 
          'Kutorginata', 'Spiriferida', 'Spiriferinida')
bins = 'StageNewOrdSplitNoriRhae20Nov2013'
gts = rev(c(lump[, 2])[[1]])
cuts = 'Chang'
bot = 'Trem'
short.data <- sort.data(bibr, payne, taxon = taxon,
                        bins = bins,
                        gts = gts,
                        cuts = cuts,
                        bot = bot)

# now just setup the data
num.samp <- nrow(short.data)

# match cohorts
con.orig <- plyr::mapvalues(short.data$orig,
                            sort(unique(as.character(short.data$orig))),
                            sort(unique(c(lump[5:(5 + 33 - 1), 4])[[1]])))
con.orig <- as.character(con.orig)
ordd <- match(con.orig, c(lump[, 4])[[1]])
con.orig <- ordd - 4
num.orig <- length(unique(con.orig))

# environmental preference
prob.epi <- qbeta(short.data$epi / (short.data$epi + short.data$off), 
                  short.data$epi.bck + 1, short.data$off.bck + 1)
env.odds <- rescale(prob.epi)
#short.data$occur <- short.data$epi + short.data$off


# gap statistic
# which are missing
to_impute <- short.data$duration <= 2
N_imp <- sum(to_impute)
N_obs <- sum(!to_impute)

gap_obs <- short.data$gap[!to_impute]

gap_obs_order <- which(!to_impute)
gap_imp_order <- which(to_impute)


#(gap_obs * (length(gap_obs) - 1) + 0.5) / length(gap_obs)

data <- list(dur = short.data$duration, 
             censored = short.data$censored,
             cohort = con.orig, 
             occupy = arm::rescale(logit(short.data$occupy)),
             env = env.odds,
             leng = arm::rescale(log(short.data$size)),
             nocc = arm::rescale(logit(short.data$rsamp)),
             to_impute = to_impute * 1,
             N_imp = N_imp,
             N_obs = N_obs,
             gap_obs = gap_obs,
             gap_obs_order = gap_obs_order,
             gap_imp_order = gap_imp_order)
data$N <- num.samp
data$O <- num.orig

with(data, {stan_rdump(list = c('N', 
                                'O',
                                'N_imp',
                                'N_obs',
                                'dur', 
                                'censored', 
                                'cohort',
                                'occupy',
                                'env', 
                                'leng',
                                'nocc',
                                'gap_obs',
                                'gap_obs_order',
                                'gap_imp_order'),
                       file = '../data/data_dump/impute_info.data.R')})
