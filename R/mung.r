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

source('../R/jade.r')
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
# argument: binning scheme being used
# argument: gts global temporal scale for the bins
# argument: cuts where the mass ext is
# arugment: bot where the top is
#
# this function cleans the data completely as document in my mung routine
# at the end, spit out the sepkoski.data file
sort.data <- function(bibr, payne, taxon = 'Rhynchonellata', 
                      bins = 'collections.stage', gts = gts,
                      cuts = 'Changhsingian',
                      bot = 'Tremadocian') {

  # bibr <- fossil
  # taxon <- 'Rhynchonellata'; bins <- 'StageNewOrdSplitNoriRhae20Nov2013'
  # gts <- rev(as.character(lump[, 2])); cuts <- 'Chang'; bot <- 'Trem'
  # i need to have good bin information, either stage 10my or fr2my
  bibr[, bins] <- as.character(bibr[, bins])
  bibr$occurrences.genus_name <- as.character(bibr$occurrences.genus_name)

  bibr <- bibr[!is.na(bibr[, bins]), ]
  bibr <- bibr[!is.na(bibr$EO_5_1_2014), ]
  bibr <- bibr[!is.na(bibr$collections.paleolngdec), ]
  bibr <- bibr[!is.na(bibr$collections.paleolatdec), ]

  bibr <- bibr[bibr$occurrences.class_name %in% taxon | 
               bibr$occurrences.order_name %in% taxon, ]
  bibr <- bibr[bibr$occurrences.genus_name %in% payne$taxon_name, ]

  # gap-iness of each stage
  dd <- split(bibr, bibr[, bins])
  dd <- llply(dd, function(x) unique(x$occurrences.genus_name))
 
  mm <- as.character(lump[, 2]) %in% names(dd)
  mm <- as.character(lump[mm, 2])
  
  dd <- dd[names(dd) %in% as.character(lump[, 2])]
  dd <- dd[names(dd)[match(mm, names(dd))]]

  gap.ratio <- c()  # gap / range through
  for(ii in seq(from = 1, to = length(dd) - 2)) {
    # before and after
    rngthr <- dd[[ii]][dd[[ii]] %in% dd[[ii + 2]]]
    # in?
    no.gap <- dd[[ii + 1]][dd[[ii + 1]] %in% rngthr]
    gap.ratio[ii] <- (length(rngthr) - length(no.gap)) / length(rngthr)
    gap.ratio[ii] <- ifelse(is.nan(gap.ratio[ii]), 0, gap.ratio[ii])
  }
  gap.ratio <- data.frame(stg = as.character(names(dd)[seq(length(dd) - 2)]), 
                          gap.ratio)
  

  #
  straight.occ <- dlply(bibr, .(occurrences.genus_name), 
                        function(x) unique(x[, bins]))
  # find out which range into the paleozoic
  too.old <- names(which(laply(straight.occ, function(x) 
                               any(which(gts %in% x) > which(gts == bot)))))
  # never in the paleozoic
  too.young <- names(which(laply(straight.occ, function(x) 
                                 all(which(gts %in% x) < which(gts == cuts)))))
  # remove those
  bibr <- bibr[!(bibr$occurrences.genus_name %in% c(too.old, too.young)), ]

  # find out which range out of the paleozoic
  straight.occ <- dlply(bibr, .(occurrences.genus_name), 
                        function(x) unique(x[, bins]))
  survivors <- names(which(laply(straight.occ, function(x) 
                                 any(which(gts %in% x) < which(gts == cuts)))))

  paleozoic <- rev(gts[which(gts == cuts):which(gts == bot)])
  bibr <- bibr[bibr[, bins] %in% paleozoic, ]

  # this section is all about finding duration
  collec.stage <- table(bibr[, bins])
  find.dur <- function(x) {
    mm <- which(gts %in% unique(x[, bins]))
    max(mm) - min(mm) + 1
  }
  # generic duration
  taxon.age <- ddply(bibr, .(occurrences.genus_name), find.dur)
  spot.fix <- split(bibr[, bins], bibr$occurrences.genus_name)
  spot.fix <- spot.fix[names(spot.fix) %in% survivors]
  spot.fix <- llply(spot.fix, function(x) 
                    max(which(gts %in% unique(x))) - which(gts %in% cuts) + 1)
  spot.fix <- melt(spot.fix)
  taxon.age[taxon.age[, 1] %in% spot.fix[, 2], 2] <- spot.fix[, 1]

  taxon.occur <- dlply(bibr, .(occurrences.genus_name), function(x) {
                       table(x[, bins])})

  rsamp <- laply(taxon.occur, function(x) sum(x) / length(x))

  taxon.gap <- c()
  for(ii in seq(length(taxon.occur))) {
    tm <- names(taxon.occur[[ii]])  # stage names
    om <- gts %in% tm  # true false list
    dm <- taxon.age[ii, 2] - 2
    taxon.gap[ii] <- (sum(om) - 2) / dm
  }
  names(taxon.gap) <- names(taxon.occur)



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

  # fixes here
  names(bibr)[names(bibr) == bins] <- 'colstage'
  ncell <- ddply(bibr, .(occurrences.genus_name, colstage), 
                 summarize, tt = length(unique(membership)))
  big.ncell <- ddply(bibr, .(colstage), summarize,
                     total = length(unique(membership)))
  
  # Chao corrected occupancy
  ncell.stage <- split(ncell, ncell[, 2])
  inc.output <- list()
  for(ii in seq(length(ncell.stage))) {
    foo <- c(sum(ncell.stage[[ii]][, 3]), ncell.stage[[ii]][, 3])
    inc.output[[ii]] <- data.frame(ncell.stage[[ii]][, 1:2], DetInc(foo))
  }
  ncell.fix <- Reduce(rbind, inc.output)


  # raw occupancy
  ncell.bygenus <- split(ncell.fix, ncell.fix$occurrences.genus_name)
  occupy <- laply(ncell.bygenus, function(x) mean(x[, 3]))
  names(occupy) <- names(ncell.bygenus)


  # want to find the number of epicontinental versus offshore
  # these are the occurrences
  onoff <- ddply(bibr, .(occurrences.genus_name), summarize,
                 epi = sum(EO_5_1_2014 == 'E'),
                 off = sum(EO_5_1_2014 == 'O'))
  # now for each stage, get the epi and off
  big.onoff <- ddply(bibr, .(colstage), summarize,
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


  # average gap percent of stages of taxon duration
  #   those before and after but not during, relative to total before after

  # get species duration along with if died in stage at/after mass extinction
  wh.stage <- llply(off, names)
  mass.ext <- cuts
  in.mass <- llply(wh.stage, function(x) x %in% mass.ext)
  censored <- laply(in.mass, function(x) {
                    o <- c()
                    if(max(which(x)) == length(x)) {
                      o <- 1
                    } else {
                      o <- 0
                    }
                    o})
  censored[which(names(off) %in% survivors)] <- 1

  orig <- laply(wh.stage, function(x) rev(gts[gts %in% x])[1])
  age.order <- llply(orig, function(x) which(gts %in% x))
  big.dead <- which(gts %in% mass.ext)
  regime <- laply(age.order, function(x) 
                  max(which(x > big.dead)))
  age.data <- cbind(taxon.age, censored, orig, regime, onoff[, -1])
  names(age.data) <- c('genus', 'duration', 'censored', 'orig', 'regime', 
                       'epi', 'off', 'epi.bck', 'off.bck')
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
  sepkoski.data$gap <- taxon.gap[match(sepkoski.data$genus, names(taxon.gap))]
  sepkoski.data$rsamp <- rsamp[match(sepkoski.data$genus, names(rsamp))]
  sepkoski.data
}


# space.time -- prep data for state-space model
# argument: bibr fossil occurrence information
# argument: taxon taxonomic group
# argument: binning scheme
# argument: gts geologic time scale
# argument: cuts where the end is
# arugment: bot where the top is
# argument: shape shape file of globe
space.time <- function(bibr, 
                       taxon = 'Rhynchonellata', 
                       bins = 'collections.stage',
                       gts = gts, 
                       cuts = 'Changhsingian',
                       bot = 'Tremadocian',
                       shape) {
  # i need to have good bin information, either stage 10my or fr2my
  bibr[, bins] <- as.character(bibr[, bins])
  bibr$occurrences.genus_name <- as.character(bibr$occurrences.genus_name)

  bibr <- bibr[!is.na(bibr[, bins]), ]
  bibr <- bibr[!is.na(bibr$collections.paleolngdec), ]
  bibr <- bibr[!is.na(bibr$collections.paleolatdec), ]
  bibr <- bibr[bibr$occurrences.class_name == taxon, ]

  paleozoic <- gts[which(gts == cuts):which(gts == bot)]
  bibr <- bibr[bibr[, bins] %in% paleozoic, ]

  # spatial analysis
  eq <- CRS("+proj=cea +lat_0=0 +lon_0=0 +lat_ts=30 +a=6371228.0 +units=m")
  globe.map <- readShapeSpatial('../data/ne_10m_coastline.shp')  # from natural earth
  proj4string(globe.map) <- eq

  spatialref <- SpatialPoints(coords = bibr[, c('collections.paleolngdec',
                                                'collections.paleolatdec')],
                              proj4string = eq)  # wgs1984.proj
  r <- raster(globe.map, nrows = 70, ncols = 34)
  sp.ras <- trim(rasterize(spatialref, r))
  membership <- cellFromXY(sp.ras, xy = bibr[, c('collections.paleolngdec', 
                                                 'collections.paleolatdec')])

  # temp vs trop
  n.tropic <- 23.5
  s.tropic <- -23.5

  north.temp <- bibr$collections.paleolatdec > n.tropic
  north.trop <- bibr$collections.paleolatdec < n.tropic &
  bibr$collections.paleolatdec > 0

  south.trop <- bibr$collections.paleolatdec > s.tropic &
  bibr$collections.paleolatdec < 0
  south.temp <- bibr$collections.paleolatdec < s.tropic

  locs <- list(north.temp, north.trop, south.trop, south.temp)

  by.loc <- llply(locs, function(x) bibr[x, ])
  occ <- llply(by.loc, function(x) {
               working <- x[, c('occurrences.genus_name', bins)]
               names(working) <- c('gen', 'stg')
               occ <- dcast(working, gen ~ stg)
               nn <- occ[, 1]
               occ <- apply(occ[, -1], 2, function(x) {
                            cc <- x > 0
                            x[cc] <- 1
                            x})
               ord <- gts %in% colnames(occ) 
               occ <- occ[, rev(gts[ord])]
               rownames(occ) <- nn
               occ})
  occ <- llply(occ, function(x) {
               if(ncol(x) < length(paleozoic)) {
                 dummy <- matrix(0, nrow = nrow(x), ncol = length(paleozoic))
                 ma <- which(paleozoic %in% colnames(x))
                 for(ii in seq(length(ma))) {
                   dummy[, ma[ii]] <- x[, ii]
                 }
                 rownames(dummy) <- rownames(x)
                 dummy
               } else {
                 x
               }})  # p/a by geologic unit for zone
  occ
}


# get the environmental occurrence for taxa at the stage level to track change
# argument: bibr data frame for all the general occurrence shit
# argument: taxonomic group
# argument: binning scheme being used
# argument: gts global temporal scale for the bins
# argument: cuts where the mass ext is
# arugment: bot where the top is
#
environ.occ <- function(bibr, 
                        taxon = 'Rhynchonellata', 
                        bins = 'StageNewOrdSplitNoriRhae20Nov2013', 
                        gts = gts, 
                        cuts = 'Chang', 
                        bot = 'Trem') {

  bibr[, bins] <- as.character(bibr[, bins])
  bibr$occurrences.genus_name <- as.character(bibr$occurrences.genus_name)

  bibr <- bibr[!is.na(bibr[, bins]), ]
  bibr <- bibr[!is.na(bibr$EO_5_1_2014), ]
  bibr <- bibr[!is.na(bibr$collections.paleolngdec), ]
  bibr <- bibr[!is.na(bibr$collections.paleolatdec), ]

  bibr <- bibr[bibr$occurrences.class_name == taxon, ]

  straight.occ <- dlply(bibr, .(occurrences.genus_name), 
                        function(x) unique(x[, bins]))
  # find out which range into the paleozoic
  too.old <- names(which(laply(straight.occ, function(x) 
                               any(which(gts %in% x) > which(gts == bot)))))
  # never in the paleozoic
  too.young <- names(which(laply(straight.occ, function(x) 
                                 all(which(gts %in% x) < which(gts == cuts)))))
  # remove those
  bibr <- bibr[!(bibr$occurrences.genus_name %in% c(too.old, too.young)), ]

  # find out which range out of the paleozoic
  straight.occ <- dlply(bibr, .(occurrences.genus_name), 
                        function(x) unique(x[, bins]))
  survivors <- names(which(laply(straight.occ, function(x) 
                                 any(which(gts %in% x) < which(gts == cuts)))))

  paleozoic <- gts[which(gts == cuts):which(gts == bot)]
  bibr <- bibr[bibr[, bins] %in% paleozoic, ]

  # now start getting to relevant occurrence information
  # split by stage
  # stage by genus
  stages <- split(bibr, bibr[, bins])
  genus <- llply(stages, function(x) split(x, x$occurrences.genus_name))
  occur <- llply(genus, function(x) llply(x, nrow))
  epicont <- llply(genus, function(x) 
                   llply(x, function(y) sum(y$EO_5_1_2014 == 'E')))

  occpat <- Map(function(a, b) cbind(total = unlist(a), epi = unlist(b)), 
                a = occur, b = epicont)
  occpat <- occpat[rev(order(match(names(occpat), gts)))]

  ambient <- laply(occpat, function(x) sum(x[, 1]))
  land <- laply(occpat, function(x) sum(x[, 2]))

  # for each observed genera
  gen <- sort(unique(bibr$occurrences.genus_name))
  oo <- list()
  for(ii in seq(length(gen))) {
    tester <- laply(occpat, function(x) any(rownames(x) == gen[ii]))
    if(sum(tester) > 1) {
      ww <- occpat[tester]
      obs <- llply(ww, function(x) x[rownames(x) == gen[ii]])
      obs.epi <- laply(obs, function(x) x[2])
      obs.off <- laply(obs, function(x) x[1] - x[2])

      back.epi <- land[tester] - obs.epi
      back.off <- ambient[tester] - land[tester] - obs.off

      oo[[ii]] <- data.frame(stg = names(obs), 
                             obs.epi, obs.off, 
                             back.epi, back.off)
    } else {
      oo[[ii]] <- NA
    }
  }
  names(oo) <- gen
  occ.filt <- Filter(function(x) length(x) > 1, oo)

  occ.filt
}
