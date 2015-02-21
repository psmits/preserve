library(reshape2)
library(plyr)
library(stringr)

source('../R/gts.r')

data.file <- list.files('../data', pattern = 'Occs')
fossil <- read.csv(paste0('../data/', data.file))

bibr <- fossil#[fossil$class_reassigned %in% c('Bivalvia', 'Rhynchonellata'), ]

# i need to have good bin information, either stage 10my or fr2my
bins <- c('collections.stage')

bibr <- bibr[!is.na(bibr[, bins]), ]

bibr <- bibr[bibr$collections.stage %in% gts, ] 
bibr$collections.stage <- as.character(bibr$collections.stage)
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
age.data <- cbind(taxon.age, censored, orig, regime)
names(age.data) <- c('genus', 'duration', 'censored', 'orig', 'regime')
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
