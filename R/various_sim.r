library(MASS)
library(arm)
library(ggplot2)
library(scales)
theme_set(theme_minimal())
theme_update(legend.position = 'none')
set.seed(420)

# basic 
draws <- mvrnorm(10000, c(0, 0), diag(c(2, 2)))
draws[, 2] <- exp(draws[, 2])
draws <- data.frame(draws)
names(draws) <- c('selectivity', 'intensity')

mvn.base <- ggplot(draws, aes(x = selectivity, y = intensity)) 
mvn.base <- mvn.base + stat_density_2d(aes(fill = ..level..), 
                                     geom = 'polygon')
mvn.base <- mvn.base + scale_fill_gradient(low = 'lightblue', 
                                         high = 'darkblue')
mvn.base <- mvn.base + scale_x_continuous(breaks = c(-1, 0, 1))
mvn.base <- mvn.base + scale_y_continuous(breaks = NULL)
mvn.base <- mvn.base + coord_cartesian(xlim = c(-3, 3))
ggsave(plot = mvn.base, filename = '../doc/talks/figure/base_hypo.png',
       width = 8, height = 6)

# range size
draws <- mvrnorm(10000, c(1, 0), diag(c(2, 2)))

draws[, 2] <- exp(draws[, 2])
draws <- data.frame(draws)
names(draws) <- c('selectivity', 'intensity')

mvn.range <- ggplot(draws, aes(x = selectivity, y = intensity)) 
mvn.range <- mvn.range + stat_density_2d(aes(fill = ..level..), 
                                         geom = 'polygon')
mvn.range <- mvn.range + scale_fill_gradient(low = 'lightblue', 
                                             high = 'darkblue')
mvn.range <- mvn.range + scale_x_continuous(breaks = c(-1, 0, 1))
mvn.range <- mvn.range + scale_y_continuous(breaks = NULL)
mvn.range <- mvn.range + coord_cartesian(xlim = c(-3, 3))
ggsave(plot = mvn.range, filename = '../doc/talks/figure/range_hypo.png',
       width = 8, height = 6)

# env preference
draws <- mvrnorm(10000, c(0, 0), diag(c(2, 2)))
draws[, 2] <- exp(draws[, 2])
draws <- data.frame(draws)
names(draws) <- c('preference', 'intensity')

mvn.env <- ggplot(draws, aes(x = preference, y = intensity)) 
mvn.env <- mvn.env + stat_density_2d(aes(fill = ..level..), 
                                     geom = 'polygon')
mvn.env <- mvn.env + scale_fill_gradient(low = 'lightblue', 
                                         high = 'darkblue')
mvn.env <- mvn.env + scale_x_continuous(breaks = c(-1, 0, 1))
mvn.env <- mvn.env + scale_y_continuous(breaks = NULL)
mvn.env <- mvn.env + coord_cartesian(xlim = c(-3, 3))
ggsave(plot = mvn.env, filename = '../doc/talks/figure/env_hypo.png',
       width = 8, height = 6)

# curavture
draws <- mvrnorm(10000, c(-1, 0), diag(c(2, 2)))
draws[, 2] <- exp(draws[, 2])
draws <- data.frame(draws)
names(draws) <- c('selectivity', 'intensity')

mvn.cur <- ggplot(draws, aes(x = selectivity, y = intensity)) 
mvn.cur <- mvn.cur + stat_density_2d(aes(fill = ..level..), 
                                     geom = 'polygon')
mvn.cur <- mvn.cur + scale_fill_gradient(low = 'lightblue', 
                                         high = 'darkblue')
mvn.cur <- mvn.cur + scale_x_continuous(breaks = c(-1, 0, 1))
mvn.cur <- mvn.cur + scale_y_continuous(breaks = NULL)
mvn.cur <- mvn.cur + coord_cartesian(xlim = c(-3, 3))
ggsave(plot = mvn.cur, filename = '../doc/talks/figure/cur_hypo.png',
       width = 8, height = 6)
