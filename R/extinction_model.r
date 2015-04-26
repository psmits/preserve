library(rstan)
library(arm)
library(parallel)
source('../R/gts.r')
source('../R/mung.r')

data.file <- list.files('../data', pattern = 'Occs')
fossil <- read.csv(paste0('../data/', data.file))
bibr <- fossil
payne <- read.table('../data/payne_bodysize/Occurrence_PaleoDB.txt',
                    header = TRUE, stringsAsFactors = FALSE)
sepkoski.data <- sort.data(bibr, payne, taxon = 'Rhynchonellata', gts)

# sepkoski.data
num.samp <- nrow(sepkoski.data)

con.orig <- match(as.character(sepkoski.data$orig), gts)
con.orig <- mapvalues(con.orig, from = unique(con.orig), 
                      to = rank(unique(con.orig)))
num.orig <- length(unique(con.orig))

con.class <- as.numeric(as.factor(sepkoski.data$class))
num.class <- length(unique(sepkoski.data$class))

# beta distribution of genus envrionmental occurrence
a <- sepkoski.data$epi
a2 <- sepkoski.data$epi.bck
b <- sepkoski.data$off
b2 <- sepkoski.data$off.bck
tax.ml <- ((a + 1) - 1) / ((a + 1) + (b + 1) - 2)
bck.ml <- ((a2 + 1) - 1) / ((a2 + 1) + (b2 + 1) - 2)
env.ml <- tax.ml - bck.ml

## beta distribution of genus lithological occurrence
#a <- sepkoski.data$car
#b <- sepkoski.data$cla
#a2 <- sepkoski.data$car.bck
#b2 <- sepkoski.data$cla.bck
#tax.ml <- ((a + 1) - 1) / ((a + 1) + (b + 1) - 2)
#bck.ml <- ((a2 + 1) - 1) / ((a2 + 1) + (b2 + 1) - 2)
#lit.ml <- tax.ml - bck.ml

# do it so i can propegate error
idv.epi <- sepkoski.data$epi
idv.off <- sepkoski.data$off
tot.epi <- sepkoski.data$epi.bck
tot.off <- sepkoski.data$off.bck
#idv.car <- sepkoski.data$car
#idv.cla <- sepkoski.data$cla
#tot.car <- sepkoski.data$car.bck
#tot.cla <- sepkoski.data$cla.bck

data <- list(duration = sepkoski.data$duration, group = con.class,
             cohort = con.orig, 
             env = env.ml,
             #lit = rescale(lit.ml),
             epi = idv.epi, off = idv.off, 
             epi.bck = tot.epi, off.bck = tot.off,
             #car = idv.car, cla = idv.cla, 
             #car.bck = tot.car, cla.bck = tot.cla,
             occupy = rescale(logit(sepkoski.data$occupy)),
             size = rescale(log(sepkoski.data$size)))

dead <- sepkoski.data$censored != 1
unc <- llply(data, function(x) x[dead])
cen <- llply(data, function(x) x[!dead])

data <- list(dur_unc = unc$duration,
             group_unc = unc$group,
             cohort_unc = unc$cohort,
             occupy_unc = unc$occupy,
             env_unc = unc$env,
             epi_unc = unc$epi,
             off_unc = unc$off,
             epi_bck_unc = unc$epi.bck,
             off_bck_unc = unc$off.bck,
#             lit_unc = unc$lit,
#             car_unc = unc$car,
#             cla_unc = unc$cla,
#             car_bck_unc = unc$car.bck,
#             cla_bck_unc = unc$cla.bck,
             size_unc = unc$size,
             N_unc = length(unc$duration),
             dur_cen = cen$duration,
             group_cen = cen$group,
             cohort_cen = cen$cohort,
             occupy_cen = cen$occupy,
             env_cen = cen$env,
             epi_cen = cen$epi,
             off_cen = cen$off,
             epi_bck_cen = cen$epi.bck,
             off_bck_cen = cen$off.bck,
#             lit_cen = cen$lit,
#             car_cen = cen$car,
#             cla_cen = cen$cla,
#             car_bck_cen = cen$car.bck,
#             cla_bck_cen = cen$cla.bck,
             size_cen = cen$size,
             N_cen = length(cen$duration))

data$samp_unc <- seq(data$N_unc)
data$samp_cen <- seq(from = data$N_unc + 1,
                     to = data$N_unc + data$N_cen,
                     by = 1)

data$N <- num.samp
data$O <- num.orig
data$C <- num.class

with(data, {stan_rdump(list = c('dur_unc', 'group_unc', 'cohort_unc', 
                                'N_unc', 'dur_cen', 'group_cen', 
                                'cohort_cen', 'N_cen', 'samp_unc', 
                                'samp_cen', 'N', 'O', 'R', 'C', 'F',
                                'occupy_unc', 'occupy_cen',
                                'env_unc', 'env_cen',
                                #'lit_unc', 'lit_cen',
                                'epi_unc', 'epi_cen', 
                                'off_unc', 'off_cen',
                                'epi_bck_unc', 'epi_bck_cen', 
                                'off_bck_unc', 'off_bck_cen', 
                                #'car_unc', 'car_cen', 
                                #'cla_unc', 'cla_cen',
                                #'car_bck_unc', 'car_bck_cen', 
                                #'cla_bck_unc', 'cla_bck_cen', 
                                'size_unc', 'size_cen'),
                       file = '../data/data_dump/fauna_info.data.R')})
