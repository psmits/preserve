library(arm)
library(parallel)
library(rjags)
library(coda)
library(rstan)
#load.module('glm')
source('../R/gts.r')
source('../R/mung.r')

data.file <- list.files('../data', pattern = 'Occs')
fossil <- read.csv(paste0('../data/', data.file))
bibr <- fossil

lump.file <- list.files('../data', pattern = 'lump')
lump <- read.csv(paste0('../data/', lump.file))

shape <- readShapeSpatial('../data/ne_10m_coastline.shp')  # from natural earth

sight <- space.time(bibr, 
                    bins = 'StageNewOrdSplitNoriRhae20Nov2013', 
                    gts = rev(as.character(lump[, 2])),
                    cuts = 'Chang',
                    bot = 'Trem',
                    shape = shape)

sizes <- laply(sight, dim)
time.bin <- sizes[1, 2]
prov.size <- sizes[, 1]

taxon <- as.numeric(as.factor(Reduce(c, llply(sight, rownames))))
ntaxa <- length(unique(taxon))

# R rows
# C columns
# T taxa
# P provinces
#
# matrix[R, C] sight;
# vector[R] taxon;
# vector[P] prov;

data <- list(R = sizes[2, 1], C = time.bin, sight = as.matrix(sight[[2]]))
jags <- with(data, {jags.model('../jags/hmm_hierarchical.jags', 
                               data = list('nyear' = C, 
                                           'nindiv' = R,
                                           'y' = sight),
                               n.chains = 4,
                               n.adapt = 1000,
                               inits = list(z = sight,
                                            p = rep(0.5, C)))})
# discard/burn-in
update(jags, 5000)
# posterior samples
post.samp <- coda.samples(jags, c('psi', 
                                  'gamma', 'phi', 'p',
                                  'z', 'turnover'),  
                          n.iter = 5000, thin = 5)
test <- coda::gelman.diag(post.samp)
