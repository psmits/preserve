library(arm)
library(plyr)
library(parallel)
library(rjags)
library(coda)
library(rstan)
load('../data/gradient_setup.rdata')
n <- 4
sight <- sight[1:n]

sizes <- laply(sight, dim)
time.bin <- as.numeric(sizes[1, 2])
prov.size <- sizes[, 1]

taxon <- as.numeric(as.factor(Reduce(c, llply(sight, rownames))))
ntaxa <- length(unique(taxon))

total.names <- sort(unique(Reduce(c, llply(sight, rownames))))

total.occur <- array(0, c(ntaxa, time.bin, length(sight)))
for(ii in seq(length(sight))) {
  total.occur[match(rownames(sight[[ii]]), total.names), , ii] <- sight[[ii]]
}

data <- list(nindiv = ntaxa, 
             nyear = time.bin, 
             nprov = length(sight), 
             y = total.occur)


save(data, file = '../data/data_dump/occurrence_data.rdata')
with(data, {dump(c('nindiv', 'nyear', 'nprov', 'y'), 
                file = '../data/data_dump/occurrence_dump.R')})

p <- array(0, dim = c(data$nyear, data$nprov))
z <- data$y
dump(c('p', 'z'), file = '../data/data_dump/hmm_inits.R')
