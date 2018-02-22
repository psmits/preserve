library(rstan)
library(rstanarm)
library(SpadeR)

library(tidyverse)

library(sp)
library(dismo)
library(raster)
library(rgdal)
library(XML)
library(maptools)

source('../R/jade.r')
#source('../R/gts.r')
#source('../R/mung.r')

# geologicl stage information
lump.file <- list.files('../data', pattern = 'lump')
lump <- read_csv(paste0('../data/', lump.file))
timeframe <- c('Trem', 'Floi', 'Dapi', 'Darr', 'Sand', 'Kati', 'Hirn', 'Ldov', 
               'Wenl', 'Ludl', 'Prid', 'Gedi', 'Sieg', 'Emsi', 'Eife', 'Give', 
               'Fras', 'Fame', 'Tour', 'Vise', 'Serp', 'Bash', 'Mosc', 'Step', 
               'Asse', 'Sakm', 'Arti', 'Kung', 'Road', 'Word', 'Capi', 'Wuchi', 
               'Chang')
timecount <- lump[lump[, 2] %in% timeframe, 1]

# raw occurrence data from pbdb
ff <- read.csv('https://paleobiodb.org/data1.2/occs/list.txt?base_name=Brachiopoda&interval=Ordovician,Permian&show=full') 

data.file <- list.files('../data', pattern = 'Occs')
fossil <- read_csv(paste0('../data/', data.file))
bibr <- fossil

# place nice
names(bibr) <- str_replace_all(names(bibr), '\\_', '\\.')
# just brachiopod orders
bo <- unique(ff$order)
bo <- bo[bo != '']

# play nice
bibr <- dplyr::mutate(bibr, stagename = StageNewOrdSplitNoriRhae20Nov2013)
bibr <- dplyr::mutate(bibr, stagenum = lump[match(stagename, lump[, 2]), 1])
bibr <- filter(bibr, !is.na(stagenum))

# survivors
bibr <- group_by(bibr, occurrences.genus.name) %>%
  dplyr::mutate(survivor = any(stagenum > max(timecount))) %>%
  ungroup()


# only species originating in timeframe
# filter down to brachiopods with info
bibr <- filter(bibr, 
               stagenum >= min(timecount), 
               stagenum <= max(timecount),
               occurrences.order.name %in% bo,  # taxa
               !is.na(collections.paleolngdec),  # location
               !is.na(collections.paleolatdec),
               !is.na(EO.5.1.2014),  # env
               stagename %in% timeframe)  # timing
bibr$stagename <- factor(bibr$stagename, levels = timeframe)



# age in # stages incl FAD, LAD
bibr <- group_by(bibr, occurrences.genus.name) %>%
  dplyr::mutate(duration = abs(max(stagenum) - min(stagenum)) + 1) %>%
  ungroup()

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
  dplyr::mutate(world = n_distinct(cell)) %>%
  group_by(occurrences.genus.name, add = TRUE) %>%
  # range size at time and more
  dplyr::mutate(range = n_distinct(cell), 
                rangeprop = range / world) %>% 
  ungroup()

# JADE correction of range proportion and more
# occurrence environment
bibr <- group_by(bibr, stagename) %>%
  dplyr::mutate(chaorange = DetInc(c(unique(world), range))) %>%
  ungroup()

# time specific occurrence types
bibr <- group_by(bibr, stagename) %>%
  dplyr::mutate(edge.world = sum(EO.5.1.2014 == 'E'),  
                off.world = sum(EO.5.1.2014 == 'O')) %>%
  # genus time specific occurrence types
  group_by(occurrences.genus.name, add = TRUE) %>%
  dplyr::mutate(edge = sum(EO.5.1.2014 == 'E'),  
                off = sum(EO.5.1.2014 == 'O')) %>%
  ungroup()



# try to move to something closer to continuous time
# connect to macrostrat
# first grab all the units collections are from
uc <- unique(bibr$collection.no)
fc <- uc[1:(length(uc) / 2)]
sc <- uc[((length(uc) / 2) + 1):length(uc)]
furl <- paste0('https://macrostrat.org/api/v2/units?cltn_id=', 
               paste0(fc, collapse = ','), 
               '&response=long&format=csv')
strat1 <- read_csv(furl)
furl <- paste0('https://macrostrat.org/api/v2/units?cltn_id=', 
               paste0(sc, collapse = ','), 
               '&response=long&format=csv')
strat2 <- read_csv(furl)
# combine into one
strat <- dplyr::union(dplyr::select(strat1, -SGp), 
                      dplyr::select(strat2, -SGp))

# then grab all the fossils from those units
furl <- paste0('https://macrostrat.org/api/v2/fossils?unit_id=', 
               paste0(unique(strat$unit_id), collapse = ','), 
               '&response=long&format=csv')
fos <- read_csv(furl)
# this gives cltn_id
fos <- dplyr::mutate(fos, collection.no = cltn_id)
ft <- left_join(bibr, fos)
ft <- dplyr::filter(ft, 
                    !is.na(t_age),
                    !is.na(b_age))

ft <- dplyr::mutate(ft, 
                    m_age = map2_dbl(t_age, b_age, ~ mean(Reduce(c, .x, .y))),
                    duration = abs(t_age - b_age))

