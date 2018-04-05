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

ss <- bibr[bibr$occurrences.genus_name %in% short.data$genus, ]
ss <- unique(ss[, c('occurrences.genus_name', 'occurrences.species_name')])
ss <- ss[order(ss[, 1]), ]
ss <- ss[!(str_detect(ss[, 2], 'sp') | str_detect(ss[, 2], '[A-Z]')), ]
ss.gen <- ddply(ss, .(occurrences.genus_name), summarize, 
                o <- length(occurrences.species_name))

# imputed model
# now just setup the data

num.samp <- nrow(short.data)

# match cohorts
con.orig <- mapvalues(short.data$orig,
                      sort(unique(as.character(short.data$orig))),
                      sort(unique(as.character(lump[5:(5 + 33 - 1), 4]))))
ordd <- match(con.orig, lump[, 4])
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
             #leng = rescale(log(short.data$size)),
             relab = rescale(log(short.data$rsamp)),
             nocc = rescale(log(short.data$occur)))


inclusion <- short.data$duration > 2

data$N <- num.samp
data$O <- num.orig

with(data, {stan_rdump(list = c('N', 
                                'O',
                                'dur', 
                                'censored', 
                                'cohort',
                                'occupy',
                                'env', 
                                'relab',
                                'nocc'),
                       file = '../data/data_dump/impute_info.data.R')})
