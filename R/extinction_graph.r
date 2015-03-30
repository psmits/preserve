library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(rstan)
library(survival)
library(stringr)
library(grid)

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 20),
             axis.title = element_text(size = 30),
             legend.text = element_text(size = 25),
             legend.title = element_text(size = 26),
             legend.key.size = unit(2, 'cm'),
             strip.text = element_text(size = 20))

fossil <- read_rdump('../data/data_dump/survival_info.data.R')
sepkoski <- read_rdump('../data/data_dump/fauna_info.data.R')

duration <- c(sepkoski$dur_unc, sepkoski$dur_cen)
condition <- c(rep(1, sepkoski$N_unc), rep(0, sepkoski$N_cen))
condition[duration == 1 & condition == 1] <- 2

emp.surv <- survfit(Surv(time = duration, time2 = duration, 
                         event = condition, type = 'interval') ~ 1)
emp.surv <- data.frame(cbind(time = emp.surv$time, surv = emp.surv$surv))
emp.surv <- rbind(c(0, 1), emp.surv)

# let's get complicated
# taxon, fauna
taxon <- c(sepkoski$group_unc, sepkoski$group_cen)
by.taxon <- split(data.frame(duration, condition), taxon)
by.taxon <- llply(by.taxon, function(x) survfit(Surv(time = x$duration, 
                                                     time2 = x$duration,
                                                     event = x$condition,
                                                     type = 'interval') ~ 1))
by.taxon <- llply(by.taxon, function(x) {
                  oo <- data.frame(time = x$time, surv = x$surv)
                  oo <- rbind(c(0, 1), oo)
                  oo})
by.taxon <- Reduce(rbind, Map(function(x, y) cbind(x, label = rep(y, nrow(x))),
                              by.taxon, names(by.taxon)))

by.taxon$fauna <- sepkoski$fauna[by.taxon$label]

taxon.surv <- ggplot(by.taxon, aes(x = time, y = surv, colour = label))
taxon.surv <- taxon.surv + geom_step(direction = 'hv')
taxon.surv <- taxon.surv + facet_wrap( ~ fauna, ncol = 1)


# cohort, regime
temporal <- c(sepkoski$cohort_unc, sepkoski$cohort_cen)
by.time <- split(data.frame(duration, condition), temporal)
by.time <- llply(by.time, function(x) survfit(Surv(time = x$duration, 
                                                   time2 = x$duration,
                                                   event = x$condition,
                                                   type = 'interval') ~ 1))
by.time <- llply(by.time, function(x) {
                 oo <- data.frame(time = x$time, surv = x$surv)
                 oo <- rbind(c(0, 1), oo)
                 oo})
by.time <- Reduce(rbind, Map(function(x, y) cbind(x, label = rep(y, nrow(x))),
                             by.time, names(by.time)))

by.time$regime <- sepkoski$regime[by.time$label]


time.surv <- ggplot(by.time, aes(x = time, y = surv, colour = label))
time.surv <- time.surv + geom_step(direction = 'hv')
time.surv <- time.surv + facet_wrap( ~ regime, ncol = 1)



# explore interactions....
inter <- interaction(taxon, temporal)
by.inter <- split(data.frame(duration, condition), inter)
by.inter <- by.inter[laply(by.inter, function(x) dim(x)[1] > 5)]
by.inter <- by.inter[laply(by.inter, function(x) any(x[, 2] == 1))]
by.inter <- llply(by.inter, function(x) survfit(Surv(time = x$duration, 
                                                   time2 = x$duration,
                                                   event = x$condition,
                                                   type = 'interval') ~ 1))
by.inter <- llply(by.inter, function(x) {
                 oo <- data.frame(time = x$time, surv = x$surv)
                 oo <- rbind(c(0, 1), oo)
                 oo})
by.inter <- Reduce(rbind, Map(function(x, y) cbind(x, label = rep(y, nrow(x))),
                             by.inter, names(by.inter)))

by.inter <- cbind(by.inter, 
                  Reduce(rbind, str_split(as.character(by.inter$label), 
                                          '\\.')))
names(by.inter) <- c('time', 'surv', 'label', 'taxon', 'cohort')

by.inter$taxon <- factor(as.numeric(as.character(by.inter$taxon)), 
                         levels = sort(unique(as.numeric(
                                  as.character(by.inter$taxon)))))
by.inter$cohort <- factor(as.numeric(as.character(by.inter$cohort)), 
                         levels = sort(unique(as.numeric(
                                  as.character(by.inter$cohort)))))
by.inter$fauna <- as.factor(sepkoski$fauna[by.inter$taxon])
by.inter$regime <- as.factor(sepkoski$regime[by.inter$cohort])

time.surv.a <- ggplot(by.inter, aes(x = time, y = surv, 
                                    group = label))
time.surv.a <- time.surv.a + geom_step(direction = 'hv', alpha = 0.5)
time.surv.a <- time.surv.a + facet_grid(regime ~ fauna)
time.surv.a <- time.surv.a + scale_colour_manual(values = cbp)

time.surv.b <- ggplot(by.inter, aes(x = time, y = surv, 
                                    group = label))
time.surv.b <- time.surv.b + geom_step(direction = 'hv', alpha = 0.5)
time.surv.b <- time.surv.b + facet_grid(fauna ~ regime)
time.surv.b <- time.surv.b + scale_colour_manual(values = cbp)
