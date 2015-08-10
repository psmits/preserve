library(arm)
library(parallel)
library(rjags)
library(coda)
library(rstan)
#load.module('glm')
source('../R/gts.r')
source('../R/mung.r')

set.seed(420)
n <- 2

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
sight <- sight[1:n]

sizes <- laply(sight, dim)
time.bin <- sizes[1, 2]
prov.size <- sizes[, 1]

taxon <- as.numeric(as.factor(Reduce(c, llply(sight, rownames))))
ntaxa <- length(unique(taxon))

total.names <- sort(unique(Reduce(c, llply(sight, rownames))))

total.occur <- array(0, c(ntaxa, time.bin, length(sight)))
for(ii in seq(length(sight))) {
  total.occur[match(rownames(sight[[ii]]), total.names), , ii] <- sight[[ii]]
}

data <- list(R = ntaxa, C = time.bin, J = length(sight), 
             sight = total.occur)
save(data, file = '../data/data_dump/occurrence_data.rdata')

p.init <- array(0, dim = c(data$C, data$J))

jags <- with(data, {jags.model('../jags/hmm_hierarchical.jags', 
                               data = list('nprov' = J,
                                           'nyear' = C, 
                                           'nindiv' = R,
                                           'y' = sight),
                               n.chains = 4,
                               n.adapt = 100,
                               inits = list(z = sight,
                                            p_norm = p.init))})
# warm-up/burn-in
update(jags, 10000)
# production
post.samp <- coda.samples(jags, c('psi', 
                                  'gamma', 'phi', 'p',
                                  'gamma_mu', 'phi_mu', 'p_mu',
                                  'gamma_sigma', 'phi_sigma', 'p_sigma',
                                  'gamma_group', 'phi_group', 'p_group',
                                  'gamma_sigma_group', 
                                  'phi_sigma_group', 
                                  'p_sigma_group',
                                  'indiv_gamma', 'indiv_phi', 'indiv_p',
                                  'indiv_gamma_sigma', 
                                  'indiv_phi_sigma', 
                                  'indiv_p_sigma',
                                  'z', 'turnover'),  
                          n.iter = 10000, thin = 10)
save(post.samp, file = '../data/mcmc_out/turnover_jags.rdata')
