library(rstan)
library(stringr)
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

con.orig <- dplyr::mapvalues(short.data$orig,
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
short.data$occur <- short.data$epi + short.data$off


data <- list(dur = short.data$duration, 
             censored = short.data$censored,
             cohort = con.orig, 
             occupy = rescale(logit(short.data$occupy)),
             env = env.odds,
             leng = rescale(log(short.data$size)),
             relab = rescale(log(short.data$rsamp)),
             nocc = rescale(log(short.data$occur)))

data$N <- num.samp
data$O <- num.orig

with(data, {stan_rdump(list = c('N', 
                                'O',
                                'dur', 
                                'censored', 
                                'cohort',
                                'occupy',
                                'env', 
                                'leng',
                                'relab',
                                'nocc'),
                       file = '../data/data_dump/impute_info.data.R')})
