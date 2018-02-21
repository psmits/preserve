library(rstan)
library(rstanarm)
library(tidyverse)
library(sp)
library(dismo)
library(raster)
library(rgdal)
library(XML)
library(maptools)
library(SpadeR)
source('../R/gts.r')
#source('../R/jade.r')
#source('../R/mung.r')

# geologicl stage information
lump.file <- list.files('../data', pattern = 'lump')
lump <- read.csv(paste0('../data/', lump.file))
timeframe <- c('Trem', 'Floi', 'Dapi', 'Darr', 'Sand', 'Kati', 'Hirn', 'Ldov', 
               'Wenl', 'Ludl', 'Prid', 'Gedi', 'Sieg', 'Emsi', 'Eife', 'Give', 
               'Fras', 'Fame', 'Tour', 'Vise', 'Serp', 'Bash', 'Mosc', 'Step', 
               'Asse', 'Sakm', 'Arti', 'Kung', 'Road', 'Word', 'Capi', 'Wuchi', 
               'Chang')

# raw occurrence data from pbdb
ff <- read.csv('https://paleobiodb.org/data1.2/occs/list.txt?base_name=Brachiopoda&interval=Ordovician,Permian&show=full', 
               stringsAsFactors = FALSE)

data.file <- list.files('../data', pattern = 'Occs')
fossil <- read.csv(paste0('../data/', data.file))
bibr <- fossil

# place nice
names(bibr) <- str_replace_all(names(bibr), '\\_', '\\.')

# just brachiopod orders
bo <- unique(ff$order)
bo <- bo[bo != '']

# play nice
bibr <- mutate(bibr, stagename = StageNewOrdSplitNoriRhae20Nov2013)
bibr <- mutate(bibr, stagenum = lump[match(stagename, lump[, 2]), 1])

# TO DO
# only species originating in timeframe



# filter down to brachiopods with info
bibr <- filter(bibr, 
               occurrences.order.name %in% bo,  # taxa
               !is.na(collections.paleolngdec),  # location
               !is.na(collections.paleolatdec),
               !is.na(EO.5.1.2014),  # env
               !is.na(stagename),
               stagename %in% timeframe)  # timing
bibr$stagename <- factor(bibr$stagename, levels = timeframe)

# put points on map
eq <- CRS("+proj=cea +lat_0=0 +lon_0=0 +lat_ts=30 +a=6371228.0 +units=m")
globe.map <- readOGR('../data/ne_10m_coastline.shp')  # from natural earth
proj4string(globe.map) <- eq
bibr <- as.data.frame(bibr)  # play nice
spatialref <- SpatialPoints(coords = bibr[, c('collections.paleolngdec', 
                                              'collections.paleolatdec')], 
                            proj4string = eq)  # wgs1984.proj
r <- raster(globe.map, nrows = 50, ncols = 50)
ras <- rasterize(spatialref, r)
# get cell # for each observation
bibr$cell <- cellFromXY(ras, 
                        xy = as.data.frame(bibr[, c('collections.paleolngdec', 
                                                    'collections.paleolatdec')]))

# how many cells are in each stage
#   world isn't sampled completely, so range size prop to world
bibr <- group_by(bibr, stagename) %>%
  dplyr::mutate(world = n_distinct(cell))

# range size at time
bibr <- group_by(bibr, stagename, occurrences.genus.name) %>%
  dplyr::mutate(range = n_distinct(cell), rangeprop = range / world) %>%
  ungroup()

# JADE correction of range proportion
bibr <- group_by(bibr, stagename) %>%
  transmute(chaorange = DetInc(c(unique(world), range))) %>%
  ungroup()

# pick up here


#ncel <- ddply(bibr, .(occurrences.genus.name, stagenum), 
#              summarize, tt = length(unique(cell)))
#ncel.stage <- split(ncel, ncel[, 2])
#sum(ncel.stage[[1]][, 3])

focus <- dplyr::select(bibr, 
                       occurrences.genus.name, 
                       age,
                       start,
                       stagename, 
                       stagenum, 
                       range, 
                       world, 
                       rangeprop)





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
