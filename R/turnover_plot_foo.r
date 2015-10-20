macro.plot <- function(posterior, 
                       time = lump,
                       label = 'Pr(z(i,t - 1) = 0 | z(i, t) = 1)',
                       filename = '../doc/gradient/figure/turnover.png',
                       subtract = FALSE,
                       process = TRUE, 
                       log = FALSE) {

  if(process) {
    est.turn <- llply(seq(data$nprov), function(y) 
                      Reduce(rbind, llply(posterior, function(x) x[[y]])))
  } else {
    est.turn <- posterior
  }
  est.mean <- melt(laply(est.turn, colMeans))
  names(est.mean) <- c('prov', 'year', 'div')
  if(subtract) {
    est.mean$div <- 1 - est.mean$div
  }

  est.turn <- Map(function(x) {
                  rownames(x) <- seq(nrow(x))
                  x}, est.turn)
  est.turn <- melt(est.turn)
  names(est.turn) <- c('sim', 'year', 'div', 'prov')
  est.turn$prov <- factor(est.turn$prov)
  if(subtract) {
    est.turn$div <- 1 - est.turn$div
  }

  # province names
  est.turn$prov <- mapvalues(est.turn$prov, unique(est.turn$prov), 
                             c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
  est.turn$prov <- factor(est.turn$prov, levels = 
                          c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
  est.mean$prov <- mapvalues(est.mean$prov, unique(est.mean$prov), 
                             c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
  est.mean$prov <- factor(est.mean$prov, levels = 
                          c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))

  # make in years
  time.slice <- time[seq(from = 6, to = data$nyear + 4), ]
  est.turn$year <- mapvalues(est.turn$year, 
                             unique(est.turn$year), time.slice[, 3])
  est.mean$year <- mapvalues(est.mean$year, 
                             unique(est.mean$year), time.slice[, 3])

  turn.est <- ggplot(est.turn, aes(x = year, y = div, group = sim))
  turn.est <- turn.est + geom_line(alpha = 0.01) 
  turn.est <- turn.est + geom_line(data = est.mean,
                                   mapping = aes(x = year, y = div, group = NULL), 
                                   colour = 'blue')
  turn.est <- turn.est + scale_x_reverse() + facet_grid(prov ~ .)
  if(log) {
    turn.est <- turn.est + scale_y_continuous(trans=log10_trans())
  }
  turn.est <- turn.est + labs(x = 'time', y = label)
  ggsave(plot = turn.est, filename = filename,
         width = 10, height = 5)
}
