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

# now just set up data
num.samp <- nrow(short.data)

con.orig <- match(as.character(short.data$orig), rev(as.character(lump[, 2])))
con.orig <- mapvalues(con.orig, from = unique(con.orig), 
                      to = rank(unique(con.orig)))
num.orig <- length(unique(con.orig))


# assemble standata
standata <- list(y = short.data$duration, 
                 status = short.data$censored,
                 cohort = con.orig)

standata$N <- num.samp
standata$O <- num.orig

# covariates
K <- 3
X <- matrix(1, ncol = K, nrow = standata$N)
X[, 2] <- rescale(logit(short.data$occupy))
X[, 3] <- rescale(log(short.data$size))
standata$X <- X
standata$K <- ncol(standata$X)

# export
with(standata, {stan_rdump(list = c('N', 
                                'O',
                                'K',
                                'y', 
                                'X',
                                'status', 
                                'cohort'
                                ),
                       file = '../data/data_dump/clean_info.data.R')})
