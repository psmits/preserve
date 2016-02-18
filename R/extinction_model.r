library(rstan)
library(arm)
library(parallel)
library(stringr)
library(ppcor)
source('../R/gts.r')
source('../R/jade.r')
source('../R/mung.r')

data.file <- list.files('../data', pattern = 'Occs')
fossil <- read.csv(paste0('../data/', data.file))
bibr <- fossil

lump.file <- list.files('../data', pattern = 'lump')
lump <- read.csv(paste0('../data/', lump.file))

payne <- read.table('../data/payne_bodysize/Occurrence_PaleoDB.txt',
                    header = TRUE, stringsAsFactors = FALSE)

taxon = c('Rhynchonellata', 'Strophomenata', 'Chileata', 'Obolellata', 
          'Kutorginata', 'Spiriferida', 'Spiriferinida')
bins = 'StageNewOrdSplitNoriRhae20Nov2013'
gts = rev(as.character(lump[, 2]))
cuts = 'Chang'
bot = 'Trem'
short.data <- sort.data(bibr, payne, taxon = taxon,
                        bins = bins,
                        gts = rev(as.character(lump[, 2])),
                        cuts = cuts,
                        bot = bot)

ss <- bibr[bibr$occurrences.genus_name %in% short.data$genus, ]
ss <- unique(ss[, c('occurrences.genus_name', 'occurrences.species_name')])
ss <- ss[order(ss[, 1]), ]
ss <- ss[!(str_detect(ss[, 2], 'sp') | str_detect(ss[, 2], '[A-Z]')), ]
ss.gen <- ddply(ss, .(occurrences.genus_name), summarize, 
                o <- length(occurrences.species_name))


# sepkoski.data
num.samp <- nrow(short.data)

con.orig <- match(as.character(short.data$orig), rev(as.character(lump[, 2])))
con.orig <- mapvalues(con.orig, from = unique(con.orig), 
                      to = rank(unique(con.orig)))
num.orig <- length(unique(con.orig))


# environmental preference
prob.epi <- qbeta(short.data$epi / (short.data$epi + short.data$off), 
                  short.data$epi.bck + 1, short.data$off.bck + 1)
env.odds <- rescale(prob.epi)



# normal data
# now just setup the data
data <- list(duration = short.data$duration, 
             cohort = con.orig, 
             occupy = rescale(logit(short.data$occupy)),
             env = env.odds,
             size = rescale(log(short.data$size)))


dead <- short.data$censored != 1
unc <- llply(data, function(x) x[dead])
cen <- llply(data, function(x) x[!dead])

data <- list(dur_unc = unc$duration,
             cohort_unc = unc$cohort,
             occupy_unc = unc$occupy,
             env_unc = unc$env,
             size_unc = unc$size,
             N_unc = length(unc$duration),
             dur_cen = cen$duration,
             cohort_cen = cen$cohort,
             occupy_cen = cen$occupy,
             env_cen = cen$env,
             size_cen = cen$size,
             N_cen = length(cen$duration))

data$N <- num.samp
data$O <- num.orig

with(data, {stan_rdump(list = c('N', 'O',
                                'N_unc', 'N_cen',
                                'dur_unc', 'cohort_unc', 
                                'dur_cen', 'cohort_cen',
                                'occupy_unc', 'occupy_cen',
                                'env_unc', 'env_cen',
                                'size_unc', 'size_cen'),
                       file = '../data/data_dump/fauna_info.data.R')})


# high graded model
# now just setup the data
high.grade <- short.data[short.data$duration > 2, ]

num.samp <- nrow(high.grade)

con.orig <- match(as.character(high.grade$orig), rev(as.character(lump[, 2])))
con.orig <- mapvalues(con.orig, from = unique(con.orig), 
                      to = rank(unique(con.orig)))
num.orig <- length(unique(con.orig))

# environmental preference
prob.epi <- qbeta(high.grade$epi / (high.grade$epi + high.grade$off), 
                  high.grade$epi.bck + 1, high.grade$off.bck + 1)
env.odds <- rescale(prob.epi)


data <- list(duration = high.grade$duration, 
             cohort = con.orig, 
             occupy = rescale(logit(high.grade$occupy)),
             samp = rescale(high.grade$gap))

dead <- high.grade$censored != 1
unc <- llply(data, function(x) x[dead])
cen <- llply(data, function(x) x[!dead])

data <- list(dur_unc = unc$duration,
             cohort_unc = unc$cohort,
             occupy_unc = unc$occupy,
             env_unc = unc$occupy,
             size_unc = unc$occupy,
             samp_unc = unc$samp,
             N_unc = length(unc$duration),
             dur_cen = cen$duration,
             cohort_cen = cen$cohort,
             occupy_cen = cen$occupy,
             env_cen = cen$occupy,
             size_cen = cen$occupy,
             samp_cen = cen$samp,
             N_cen = length(cen$duration))

data$N <- num.samp
data$O <- num.orig
data$T <- 3

data$cohort_cen <- data$cohort_cen[!(data$dur_cen == 3)]
data$occupy_cen <- data$occupy_cen[!(data$dur_cen == 3)]
data$samp_cen <- data$samp_cen[!(data$dur_cen == 3)]
data$env_cen <- data$env_cen[!(data$dur_cen == 3)]
data$size_cen <- data$size_cen[!(data$dur_cen == 3)]
data$dur_cen <- data$dur_cen[!(data$dur_cen == 3)]
data$N_cen <- length(data$dur_cen)

with(data, {stan_rdump(list = c('N', 'O', 'T',
                                'N_unc', 'N_cen',
                                'dur_unc', 'cohort_unc', 
                                'dur_cen', 'cohort_cen',
                                'occupy_unc', 'occupy_cen',
                                'env_unc', 'env_cen',
                                'size_unc', 'size_cen',
                                'samp_unc', 'samp_cen'),
                       file = '../data/data_dump/high_info.data.R')})



# imputed model
# now just setup the data

num.samp <- nrow(short.data)

con.orig <- match(as.character(short.data$orig), rev(as.character(lump[, 2])))
con.orig <- mapvalues(con.orig, from = unique(con.orig), 
                      to = rank(unique(con.orig)))
num.orig <- length(unique(con.orig))

# environmental preference
prob.epi <- qbeta(short.data$epi / (short.data$epi + short.data$off), 
                  short.data$epi.bck + 1, short.data$off.bck + 1)
env.odds <- rescale(prob.epi)


data <- list(dur = short.data$duration, 
             censored = short.data$censored,
             cohort = con.orig, 
             occupy = rescale(logit(short.data$occupy)),
             env = env.odds,
             leng = rescale(log(short.data$size)))


inclusion <- (!(is.nan(short.data$gap))) & (!(is.infinite(short.data$gap)))

data$samp_obs <- short.data$gap[inclusion]
data$N <- num.samp
data$O <- num.orig
data$N_obs <- length(data$samp_obs)
data$N_imp <- data$N - data$N_obs
data$inclusion <- inclusion * 1


with(data, {stan_rdump(list = c('N', 
                                'O',
                                'N_obs', 'N_imp',
                                'censored', 
                                'inclusion',
                                'dur', 
                                'cohort',
                                'occupy',
                                'env', 
                                'leng', 
                                'samp_obs'),
                       file = '../data/data_dump/impute_info.data.R')})
