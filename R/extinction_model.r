library(rstan)
library(arm)
library(parallel)
library(stringr)
library(ppcor)
source('../R/gts.r')
source('../R/mung.r')

data.file <- list.files('../data', pattern = 'Occs')
fossil <- read.csv(paste0('../data/', data.file))
bibr <- fossil

lump.file <- list.files('../data', pattern = 'lump')
lump <- read.csv(paste0('../data/', lump.file))

payne <- read.table('../data/payne_bodysize/Occurrence_PaleoDB.txt',
                    header = TRUE, stringsAsFactors = FALSE)
short.data <- sort.data(bibr, payne, taxon = 'Rhynchonellata', 
                        bins = 'StageNewOrdSplitNoriRhae20Nov2013', 
                        gts = rev(as.character(lump[, 2])),
                        cuts = 'Chang',
                        bot = 'Trem')

ss <- bibr[bibr$occurrences.genus_name %in% short.data$genus, ]
ss <- unique(ss[, c('occurrences.genus_name', 'occurrences.species_name')])
ss <- ss[order(ss[, 1]), ]
ss <- ss[!(str_detect(ss[, 2], 'sp') | str_detect(ss[, 2], '[A-Z]')), ]
ss.gen <- ddply(ss, .(occurrences.genus_name), summarize, 
                o <- length(occurrences.species_name))

#tt <- data.frame(dur = short.data$duration, 
#                 obs = short.data$occupy,
#                 avg = short.data$epi + short.data$off / short.data$duration)
#cor(tt$dur, tt$obs)
#cor(tt$dur, tt$avg)
#cor(tt$obs, tt$avg)
#pcor(tt)


# sepkoski.data
num.samp <- nrow(short.data)

con.orig <- match(as.character(short.data$orig), rev(as.character(lump[, 2])))
con.orig <- mapvalues(con.orig, from = unique(con.orig), 
                      to = rank(unique(con.orig)))
num.orig <- length(unique(con.orig))


# do it so i can propegate error
idv.epi <- short.data$epi
idv.off <- short.data$off

data <- list(duration = short.data$duration, 
             cohort = con.orig, 
             epi = idv.epi, 
             off = idv.off, 
             occupy = rescale(logit(short.data$occupy)),
             size = rescale(log(short.data$size)))

dead <- short.data$censored != 1
unc <- llply(data, function(x) x[dead])
cen <- llply(data, function(x) x[!dead])

data <- list(dur_unc = unc$duration,
             cohort_unc = unc$cohort,
             occupy_unc = unc$occupy,
             epi_unc = unc$epi,
             off_unc = unc$off,
             size_unc = unc$size,
             N_unc = length(unc$duration),
             dur_cen = cen$duration,
             cohort_cen = cen$cohort,
             occupy_cen = cen$occupy,
             epi_cen = cen$epi,
             off_cen = cen$off,
             size_cen = cen$size,
             N_cen = length(cen$duration))

data$N <- num.samp
data$O <- num.orig

with(data, {stan_rdump(list = c('N', 'O',
                                'N_unc', 'N_cen',
                                'dur_unc', 'cohort_unc', 
                                'dur_cen', 'cohort_cen',
                                'occupy_unc', 'occupy_cen',
                                'epi_unc', 'epi_cen', 
                                'off_unc', 'off_cen',
                                'size_unc', 'size_cen'),
                       file = '../data/data_dump/fauna_info.data.R')})
