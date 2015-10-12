library(plyr)
library(coda)
library(arm)
library(stringr)
library(parallel)
source('../R/turnover_functions.r')

cod <- mclapply(1:4, function(x){
                oo <- paste0('../data/mcmc_out/CODAchain', x, '.txt')
                cc <- read.coda(output.file = oo, 
                                index.file = '../data/mcmc_out/CODAindex.txt')
                cc}, mc.cores = detectCores())
post.samp <- mcmc.list(cod)

# read in 4 coda files
# make multi mcmc object
# process into stan like format because easier to read
#load('../data/data_dump/occurrence_data.rdata')  # data
load('../data/mcmc_out/turnover_jags.rdata')  # post.samp

post <- process.coda(post.samp, data)
save(post, file = '../data/mcmc_out/turnover_custom.rdata')
