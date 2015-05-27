library(reshape2)
library(arm)
library(plyr)
library(stringr)
library(raster)
library(igraph)
library(sp)
library(dismo)
library(XML)
library(maptools)
library(foreign)
library(rgdal)

source('../R/clean_funcs.r')
source('../R/gts.r')
sepkoski <- list(cambrian = c('Trilobita', 'Polychaeta', 'Tergomya', 
                              'Lingulata'),
                 paleozoic = c('Rhynchonellata', 'Crinoidea', 'Ostracoda', 
                               'Cephalopoda', 'Anthozoa', 'Cyclocystoidea', 
                               'Asteroidea', 'Ophiuroidea'),
                 modern = c('Gastropoda', 'Bivalvia', 'Osteichtyes', 
                            'Malacostraca', 'Echinoidea', 'Gymnolaemata', 
                            'Demospongea', 'Chondrichthyes'))

# sort.data -- prep data for survival analysis
# start function here, runs all the way to the end of the file
# argument: bibr data frame for all the general occurrence shit
# argument: payne body size data
# argument: taxonomic group
# argument: gts global temporal scale by stages
# (argument: temporal unit; currently it is set up for pre-P/T boundary)
#
# this function cleans the data completely as document in my mung routine
# at the end, spit out the sepkoski.data file
sort.data <- function(bibr, payne, taxon = 'Rhynchonellata', gts = gts) {

  # i need to have good bin information, either stage 10my or fr2my
  bins <- c('collections.stage')
  bibr$collections.stage <- as.character(bibr$collections.stage)
  bibr$occurrences.genus_name <- as.character(bibr$occurrences.genus_name)
  
  bibr <- bibr[!is.na(bibr[, bins]), ]
  bibr <- bibr[!is.na(bibr$EO_5_1_2014), ]
  bibr <- bibr[!is.na(bibr$collections.paleolngdec), ]
  bibr <- bibr[!is.na(bibr$collections.paleolatdec), ]
  bibr <- bibr[bibr$collections.stage %in% gts, ] 

  paleozoic <- gts[which(gts == 'Changhsingian'):length(gts)]
  bibr <- bibr[bibr$collections.stage %in% paleozoic, ]

  bibr <- bibr[bibr$occurrences.class_name == taxon, ]
  bibr <- bibr[bibr$occurrences.genus_name %in% payne$taxon_name, ]

  # lithology
  # high chance of removing occurrences
  bibr$collections.lithology1 <- as.character(bibr$collections.lithology1)
  bibr <- bibr[!is.na(bibr$collections.lithology1), ]
  lith <- bibr$collections.lithology1
  lith <- gsub(pattern = '[\\"?]', replacement = '', lith, perl = TRUE)
  bibr$collections.lithology1 <- clean.lith(lith)
  bibr <- bibr[!(bibr$collections.lithology1 %in% c('', 'mixed')), ]

  # this section is all about finding duration
  collec.stage <- table(bibr$collections.stage)
  find.dur <- function(x) {
    mm <- which(gts %in% unique(x$collections.stage))
    max(mm) - min(mm) + 1
  }
  # generic duration
  taxon.age <- ddply(bibr, .(occurrences.genus_name), find.dur)
  taxon.occur <- dlply(bibr, .(occurrences.genus_name), function(x) {
                       table(x$collections.stage)})
  taxon.nstage <- unlist(llply(taxon.occur, length))

  nzero <- taxon.age[, 2] - taxon.nstage

  # this is about geographic range size
  eq <- CRS("+proj=cea +lat_0=0 +lon_0=0 +lat_ts=30 +a=6371228.0 +units=m")
  globe.map <- readShapeSpatial('../data/ne_10m_coastline.shp')  # from natural earth
  proj4string(globe.map) <- eq
  spatialref <- SpatialPoints(coords = bibr[, c('collections.lngdec', 
                                                'collections.latdec')], 
                              proj4string = eq)  # wgs1984.proj
  r <- raster(globe.map, nrows = 70, ncols = 34)
  sp.ras <- rasterize(spatialref, r)
  bibr$membership <- cellFromXY(sp.ras, xy = bibr[, c('collections.lngdec', 
                                                      'collections.latdec')])
  ncell <- ddply(bibr, .(occurrences.genus_name, collections.stage), 
                 summarize, tt = length(unique(membership)))
  big.ncell <- ddply(bibr, .(collections.stage), summarize,
                     total = length(unique(membership)))

  ncell.bygenus <- split(ncell, ncell$occurrences.genus_name)
  occupy <- llply(ncell.bygenus, function(x) {
                  xx <- match(x[, 2], big.ncell$collections.stage)
                  xx <- x$tt / big.ncell[xx, 2]
                  mean(xx)})


  # want to find the number of epicontinental versus offshore
  # these are the occurrences
  onoff <- ddply(bibr, .(occurrences.genus_name), summarize,
                 epi = sum(EO_5_1_2014 == 'E'),
                 off = sum(EO_5_1_2014 == 'O'))
  # now for each stage, get the epi and off
  big.onoff <- ddply(bibr, .(collections.stage), summarize,
                     epi = sum(EO_5_1_2014 == 'E'),
                     off = sum(EO_5_1_2014 == 'O'))
  for(ii in seq(length(taxon.occur))) {
    app <- which(gts %in% names(taxon.occur[[ii]]))
    wh <- gts[seq(min(app), max(app))]
    background <- big.onoff[big.onoff[, 1] %in% wh, ]
    epi.back <- sum(background$epi) - onoff[ii, 2]
    off.back <- sum(background$off) - onoff[ii, 3]
    onoff[ii, 4] <- epi.back
    onoff[ii, 5] <- off.back
  }

  # do the above based on lithology
  litho <- ddply(bibr, .(occurrences.genus_name), summarize,
                 carbonate = sum(collections.lithology1 == 'carbonate'),
                 clastic = sum(collections.lithology1 == 'clastic'))
  # now for each stage, get the epi and off
  big.litho <- ddply(bibr, .(collections.stage), summarize,
                     carbonate = sum(collections.lithology1 == 'carbonate'),
                     clastic = sum(collections.lithology1 == 'clastic'))
  for(ii in seq(length(taxon.occur))) {
    app <- which(gts %in% names(taxon.occur[[ii]]))
    wh <- gts[seq(min(app), max(app))]
    background <- big.litho[big.litho[, 1] %in% wh, ]
    car.back <- sum(background$car) - litho[ii, 2]
    cla.back <- sum(background$cla) - litho[ii, 3]
    litho[ii, 4] <- car.back
    litho[ii, 5] <- cla.back
  }


  # the number of collections to offset each observation by 
  off <- list()
  for(ii in seq(length(taxon.occur))) {
    app <- which(gts %in% names(taxon.occur[[ii]]))
    wh <- gts[seq(min(app), max(app))]
    off[[ii]] <- collec.stage[names(collec.stage) %in% wh]
  }
  names(off) <- names(taxon.occur)

  # occurrences for each stage, inclusive, of generic duration
  for(ii in seq(length(taxon.occur))) {
    blank <- names(off[[ii]])
    keep <- rep(0, length(blank))
    keep[which(blank %in% names(taxon.occur[[ii]]))] <- taxon.occur[[ii]]
    taxon.occur[[ii]] <- keep
  }
  occurs <- lapply(taxon.occur, length)

  # get species duration along with if died in stage before/of mass extinction
  wh.stage <- llply(off, names)
  mass.ext <- c('Masstrichtian', 'Rhaetian', 'Changhsingian', 
                'Frasnian', 'Hirnantian', 'Calabrian')
  in.mass <- llply(wh.stage, function(x) x %in% mass.ext)
  censored <- laply(in.mass, function(x) {
                    o <- c()
                    if(max(which(x)) == length(x)) {
                      o <- 1
                    } else {
                      o <- 0
                    }
                    o})
  orig <- laply(wh.stage, function(x) rev(gts[gts %in% x])[1])
  age.order <- llply(orig, function(x) which(gts %in% x))
  big.dead <- which(gts %in% mass.ext)
  regime <- laply(age.order, function(x) 
                  max(which(x > big.dead)))
  age.data <- cbind(taxon.age, censored, orig, regime, onoff[, -1], litho[, -1])
  names(age.data) <- c('genus', 'duration', 'censored', 'orig', 'regime', 
                       'epi', 'off', 'epi.bck', 'off.bck',
                       'car', 'cla', 'car.bck', 'cla.bck')
  # need to retain the class stuff too

  # split based on class
  ords <- unique(bibr[, c('class_reassigned', 'occurrences.genus_name')])
  ords <- ords[!is.na(ords[, 1]), ]
  ords <- apply(ords, 2, as.character)
  ords <- data.frame(ords)
  class.mem <- split(ords, ords[, 1])
  class.mem <- llply(class.mem, function(x) {
                     names(taxon.occur) %in% x[, 2]})

  ords <- ords[order(ords[, 1]), ]
  ords[, 1] <- as.character(ords[, 1])

  occ.gen <- melt(taxon.occur)
  occ.gen <- occ.gen[occ.gen[, 2] %in% ords[, 2], ]

  off.melt <- melt(off)
  off.melt <- off.melt[off.melt[, 3] %in% ords[, 2], ]

  ords <- ords[match(occ.gen[, 2], ords[, 2]), ]

  dat.full <- cbind(occ.gen, ords[, 1], offset = off.melt$value)
  names(dat.full) <- c('count', 'genus', 'order', 'offset')

  # finish up the duration stuff
  age.data <- age.data[age.data$genus %in% ords[, 2], ]
  age.data$class <- ords[match(age.data$genus, ords[, 2]), 1]

  # get the subset that corresponds to the sepkoski fauna
  fauna <- age.data$class
  ws <- llply(sepkoski, function(x) age.data$class %in% x)

  for(ii in seq(length(ws))) {
    fauna[ws[[ii]]] <- names(sepkoski)[ii]
  }
  sepkoski.data <- cbind(age.data, fauna)[fauna %in% names(sepkoski), ]
  sepkoski.data$fauna <- as.character(sepkoski.data$fauna)
  sepkoski.data$occupy <- unlist(occupy[match(sepkoski.data$genus, names(occupy))])
  sepkoski.data$size <- payne$size[match(sepkoski.data$genus, payne$taxon_name)]
  sepkoski.data
}


# space.time -- prep data for state-space model
space.time <- function(bibr, taxon = 'Rhynchonellata', gts = gts, shape) {
  # i need to have good bin information, either stage 10my or fr2my
  #taxon <- 'Rhynchonellata'
  #bibr <- fossil
  bins <- c('collections.stage')
  bibr$collections.stage <- as.character(bibr$collections.stage)
  bibr$occurrences.genus_name <- as.character(bibr$occurrences.genus_name)
  bibr <- bibr[bibr$occurrences.class_name == taxon, ]
 
  to.rm <- apply(is.na(bibr[, c('occurrences.genus_name', 
                                'collections.stage', 
                                'collections.paleolngdec', 
                                'collections.paleolatdec')]), 
                 1, any)
  bibr <- bibr[!to.rm, ]
  bibr <- bibr[bibr$collections.stage %in% gts, ] 
  paleozoic <- gts[which(gts == 'Changhsingian'):length(gts)]
                   #which(gts == 'Asselian')]
  bibr <- bibr[bibr$collections.stage %in% paleozoic, ]

  # to do this by province 
  # geographic position
  eq <- CRS("+proj=cea +lat_0=0 +lon_0=0 +lat_ts=30 +a=6371228.0 +units=m")
  proj4string(shape) <- eq
  spatialref <- SpatialPoints(coords = bibr[, c('collections.paleolngdec', 
                                                'collections.paleolatdec')], 
                              proj4string = eq)  # wgs1984.proj
  r <- raster(shape, nrows = 50, ncols = 50)
  sp.ras <- rasterize(spatialref, r)
  bibr$membership <- cellFromXY(sp.ras, xy = bibr[, c('collections.lngdec', 
                                                      'collections.latdec')])
  # identify bioprovinces
  cooc <- bibr[, c('occurrences.genus_name', 'membership')]
  g <- graph.data.frame(cooc, directed = FALSE)
  V(g)$type <- V(g)$name %in% unique(cooc[, 2])
  mapcom <- infomap.community(g, nb.trials = 1000)
  
  cg <- contract.vertices(g, membership(mapcom))
  E(cg)$weight <- 1
  cg <- simplify(cg, remove.loops = FALSE)
  cgcom <- infomap.community(cg, nb.trials = 1000)
  
  # what community do the originals belong to?
  # what community do the second level belong to?
  members <- taxon.member <- loc.member <- list()
  taxa <- which(str_detect(V(g)$name, '[A-Za-z]'))
  loc <- which(str_detect(V(g)$name, '[0-9]'))
  for(ii in seq(length(communities(cgcom)))) {
    members[[ii]] <- unlist(communities(mapcom)[communities(cgcom)[[ii]]])
    taxon.member[[ii]] <- members[[ii]][members[[ii]] %in% taxa]
    loc.member[[ii]] <- members[[ii]][members[[ii]] %in% loc]
  }
  # assign taxa to provinces
  second.member <- rep(seq(length(members)), times = laply(members, length))
  cg2 <- contract.vertices(cg, membership(cgcom))
  E(cg2)$weight <- 1
  cg3 <- simplify(cg2, remove.loops = TRUE)
  E(cg3)$weight <- degree(cg2)
  # TODO

  # each row corresponds to a taxons occurrence in a province
  # have province indicator
  # have taxon indicator

  working <- bibr[, c('occurrences.genus_name', 'collections.stage')]
  occ <- dcast(working, occurrences.genus_name ~ collections.stage)
  occ <- apply(occ[, -1], 2, function(x) {
               cc <- x > 0
               x[cc] <- 1
               x})
  ord <- match(colnames(occ), paleozoic)
  ord <- mapvalues(ord, from = unique(ord), to = rank(unique(ord)))
  occ <- occ[, ord]
  occ
}
