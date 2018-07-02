library(rstan)
library(stringr)
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

#ss <- bibr[bibr$occurrences.genus_name %in% short.data$genus, ]
#ss <- unique(ss[, c('occurrences.genus_name', 'occurrences.species_name')])
#ss <- ss[order(ss[, 1]), ]
#ss <- ss[!(str_detect(ss[, 2], 'sp') | str_detect(ss[, 2], '[A-Z]')), ]
#ss.gen <- ddply(ss, .(occurrences.genus_name), summarize, 
#                o <- length(occurrences.species_name))

# imputed model
# now just setup the data
# match cohorts

con.orig <- plyr::mapvalues(short.data$orig,
                            sort(unique(as.character(short.data$orig))),
                            sort(as.character(lump[5:(5 + 33 - 1), 4])))
con.orig <- as.character(con.orig)
ordd <- match(con.orig, as.character(lump[, 4]))
con.orig <- ordd - 4
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
             leng = rescale(log(short.data$size)),
             relab = rescale(log(short.data$rsamp)))


inclusion <- short.data$duration > 2

data$samp_obs <- short.data$gap[inclusion]

data$obs_ord <- which(inclusion)
data$imp_ord <- which(!inclusion)


data$N <- length(data$dur)
data$O <- num.orig
data$N_obs <- length(data$samp_obs)
data$N_imp <- data$N - data$N_obs
data$inclusion <- inclusion * 1

# suggested by betareg manual, which has the citation
num <- data$samp_obs * (length(data$samp_obs) - 1) + 0.5
data$samp_obs <- num / length(data$samp_obs)

data$K <- 5
data$L <- 1

# export for stan
with(data, {stan_rdump(list = c('N', 
                                'O',
                                'K',
                                'L',
                                'N_obs', 'N_imp',
                                'dur', 
                                'censored', 
                                'inclusion',
                                'obs_ord',
                                'imp_ord',
                                'cohort',
                                'occupy',
                                'env', 
                                'leng',
                                'relab',
                                'samp_obs'),
                       file = '../data/data_dump/impute_info.data.R')})
