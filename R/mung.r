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

for(ii in seq(length(taxon.occur))) {
  blank <- names(off[[ii]])
  keep <- rep(0, length(blank))
  keep[which(blank %in% names(taxon.occur[[ii]]))] <- taxon.occur[[ii]]
  taxon.occur[[ii]] <- keep
}
occurs <- lapply(taxon.occur, length)


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
