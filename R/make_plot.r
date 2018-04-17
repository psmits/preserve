# start the whole plotting sequence
#########################
source('../R/multiplot.r')
source('../R/borrow_plotcorr.r')
source('../R/plot_foo.r')
posterior.plots <- function(data, 
                            wei.fit, 
                            npred, 
                            lump, 
                            name = 'cweib', 
                            left = FALSE) {
  # data setup
  # HERE
  coh <- c(data$cohort)
  rage <- c(data$occupy)
  envs <- c(data$env)
  size <- c(data$size)
  duration <- c(data$dur)

  cn <- lump[5:(5 + 33 - 1), 4] # time unit names
  cr <- purrr::map2_chr(cn, seq(length(cn)), ~ paste0(.y, '. ', .x)) # give number
  mvgr <- mapvalues(data$cohort, sort(unique(data$cohort)), cr)
  mvgr <- factor(mvgr, levels = cr) 

  # basic posterior predictive checks
  #   points based
  #   distribution based
  plot_ppcbasics(wei.fit = wei.fit, 
                 data = data, 
                 lump = lump, 
                 mvgr = mvgr, 
                 cn = cn, 
                 name = name)


  # data distributions would be cool, specifically in the case of the imputation

  ss <- sample(length(wei.fit$lp__), 12)

  # lets make survival curves
  # non-parametric estimate of the survival curve just given durations
  #   some durations are right-censored
  condition <- (data$censored == 0) * 1
  duration1 <- ifelse(condition == 1 & duration == 1, NA, duration)
  duration2 <- ifelse(condition == 0, NA, duration)

  if(!left) {
    emp.surv <- survfit(Surv(time = duration, 
                             event = condition, 
                             type = 'right') ~ 1)
  } else if(left) {
    emp.surv <- survfit(Surv(time = duration1, 
                             time2 = duration2,
                             type = 'interval2') ~ 1)
  }
  emp.surv <- data.frame(time = emp.surv$time, surv = emp.surv$surv, 
                         lower = emp.surv$lower, upper = emp.surv$upper)

  # i have two ways of estimating the survival functions from posterior
  #   actual S(t) estimates
  #   estimated durations under discrete Weibull
  wr <- wei.fit$y_tilde[sample(nrow(wei.fit$y_tilde), 100), ]
  if(!left) {
    wei.surv <- apply(wr, 1, function(x) 
                      survfit(Surv(x, event = condition) ~ 1))
  } else if(left) {
    wei.surv <- apply(wr, 1, function(x)  {
                        duration1 <- ifelse(condition == 1 & duration == 1, 
                                            NA, x)
                        duration2 <- ifelse(condition == 0, NA, x)
                        survfit(Surv(time = duration1, 
                                     time2 = duration2, 
                                     type = 'interval2') ~ 1)})
  }

  wei.surv <- llply(wei.surv, function(x) {
                      y <- data.frame(time = x$time, surv = x$surv)
                      y})
  wei.surv <- Reduce(rbind, Map(function(x, y) {
                                  x$group <- y
                                  x}, 
                                  x = wei.surv, 
                                  y = seq(length(wei.surv))))


  # naming/legacy code issue
  sim.surv <- wei.surv

  # fit model
  surv.plot <- ggplot(emp.surv, aes(x = time, y = surv))
  surv.plot <- surv.plot + geom_line(data = sim.surv, 
                                     aes(x = time, y = surv, group = group),
                                     colour = 'black', alpha = 0.05)
  surv.plot <- surv.plot + geom_line(colour = 'blue', size = 1.2)
  surv.plot <- surv.plot + geom_ribbon(mapping = aes(ymin = lower, ymax = upper), 
                                       fill = 'blue', size = 1.2, alpha = 0.3)
  surv.plot <- surv.plot + coord_cartesian(xlim = c(-0.5, max(duration)))
  surv.plot <- surv.plot + labs(x = 'Duration (t)', 
                                y = 'Pr(t < T)')
  surv.plot <- surv.plot + theme(axis.title = element_text(size = 25))
  #surv.plot <- surv.plot + scale_y_continuous(trans=log10_trans(),
  #                                            breaks = c(0.01, 0.1, 0.5, 1))
  ggsave(surv.plot, filename = paste0('../doc/figure/survival_curves_', 
                                      name, '.pdf'),
         width = 6, height = 5, dpi = 600)


  # survival curves by grouping
  #duration
  #condition
  #coh
  survdat <- data.frame(duration, condition, coh, duration1, duration2)
  survdat_g <- split(survdat, survdat$coh)
  if(!left) {
    survfit_g <- purrr::map(survdat_g, ~ survfit(Surv(time = .x$duration, 
                                                      event = .x$condition, 
                                                      type = 'right') ~ 1))
  } else if(left) {
    survfit_g <- purrr::map(survdat_g, ~ survfit(Surv(time = .x$duration1, 
                                                      time2 = .x$duration2,
                                                      type = 'interval2') ~ 1))
  }
  survfit_g <- purrr::imap(survfit_g, ~ data.frame(time = .x$time,
                                                   surv = .x$surv, 
                                                   group = .y))
#  survfit_g <- purrr::map2(survfit_g, names(survfit_g), function(x, y) {
#                             y <- data.frame(time = x$time,
#                                             surv = x$surv, 
#                                             group = y)
#                             y})
  sfg <- purrr::reduce(survfit_g, rbind)

  # now for the simulations
  gg <- sample(nrow(wei.fit$y_tilde), 100)
  wr <- wei.fit$y_tilde[gg, ]
  wrl <- as.list(as.data.frame(t(wr)))
  if(!left) {
    wrl <- purrr::map(wrl, ~ data.frame(time = .x, condition))
    wrl_g <- purrr::map(wrl, function(x) {
                          a <- split(x, coh)
                          purrr::map(a, ~ survfit(Surv(time = .x$time,
                                                       event = .x$condition)
                          ~ 1))})
  } else if(left) {
    wrl <- purrr::map(wrl, function(x) {
                        dr1 <- ifelse(duration == 1 & condition == 1, NA, x)
                        dr2 <- ifelse(condition == 0, NA, x)
                        data.frame(time = dr1, time2 = dr2)})
    x <- wrl[[1]]
    wrl_g <- purrr::map(wrl, function(x) {
                          a <- split(x, coh)
                          purrr::map(a, ~ survfit(Surv(time = .x$time,
                                                       time2 = .x$time2,
                                                       type = 'interval2') 
                          ~ 1))})
  }
  wrl_f <- purrr::map(wrl_g, function(x) {
                        purrr::imap(x, ~ data.frame(time = .x$time, 
                                                    surv = .x$surv,
                                                    group = .y))
                        #purrr::map2(x, names(x), function(a, b) {
                        #              y <- data.frame(time = a$time, 
                        #                              surv = a$surv,
                        #                              group = b)
                        #              y})})
  wrl_f <- purrr::map(wrl_f, ~ purrr::reduce(.x, rbind))
  wrl_f <- bind_rows(wrl_f, .id = 'sim')
  wrl_f$group <- as.character(wrl_f$group)


  mvg <- levels(mvgr)
  unique(sfg$group)
  sfg$group <- plyr::mapvalues(sfg$group,
                               from = unique(sfg$group),
                               to = mvg)
  sfg$group <- factor(sfg$group, levels = mvg)
  wrl_f$group <- plyr::mapvalues(wrl_f$group, 
                                 from = unique(wrl_f$group),
                                 to = mvg)
  wrl_f$group <- factor(wrl_f$group, levels = mvg)

  # facet-d by group
  sgg <- ggplot(sfg, aes(x = time, y = surv))
  sgg <- sgg + geom_line(data = wrl_f, mapping = aes(x = time, y = surv, group = sim), 
                         alpha = 0.1)
  sgg <- sgg + geom_line(colour = 'blue')
  sgg <- sgg + facet_wrap(~ group, strip.position = 'bottom')
  sgg <- sgg + coord_cartesian(xlim = c(0, 30))
  sgg <- sgg + theme(strip.text = element_text(size = 10))
  sgg <- sgg + labs(x = 'Duration (geological stages)', 
                    y = 'P(T > t)')
  ggsave(filename = '../doc/figure/ppc_surv_coh.png', sgg,
         width = 10.5, height = 8, dpi = 600)

  sgg <- sgg + scale_y_continuous(trans = log_trans())
  ggsave(filename = '../doc/figure/ppc_surv_coh_log.png', sgg,
         width = 10.5, height = 8, dpi = 600)



  wei.covcor <- get.covcor(stanfit = wei.fit, npred = npred)

  # make a plot of the correlation matrix
  col1 <- colorRampPalette(c("red", "white", "blue"))
  col1<- col1(200)
  col2 <- colorRampPalette(c("grey", "white", "grey"))
  col2 <- col2(200)

  png(file = paste0('../doc/figure/wei_cor_heatmap_', name, '.png'), 
      width = 1500, height = 1500)
  my.plotcorr(wei.covcor[[1]], upper.panel = 'number', 
              col = col1[((wei.covcor[[1]] + 1)/2) * 200], 
              mar = rep(0, 4), cex = 4, cex.lab = 4.5,
              npred = npred)
  dev.off()

  png(file = paste0('../doc/figure/wei_cor_bw_', name, '.png'), 
      width = 1500, height = 1500)
  my.plotcorr(wei.covcor[[1]], upper.panel = 'number', 
              col = col2[((wei.covcor[[1]] + 1)/2) * 200], 
              mar = rep(0, 4), cex = 4, cex.lab = 4.5,
              npred = npred)
  dev.off()


  # plots of the covariate effects estimates over time
  plot_coveffect(wei.fit = wei.fit,
                 npred = npred, 
                 data = data,
                 lump = lump,
                 name = name)



  rang <- c(data$occupy)
  if(npred == 6) {
    val2 <- quantile(rang, c(0.2, 0.5, 0.8)) 
    # very important variable for all plots that follow!
    # need to loop through each value and make a diff plot for each
    type <- c('low', 'med', 'high')
  } else if(npred == 5) {
    val2 <- quantile(rang, 0.5)
    # very important variable for all plots that follow!
    # need to loop through each value and make a diff plot for each
    type <- 'med'
  }

  sam <- sample(length(wei.fit$lp__))

  # plots of the effects of environmental preference 
  #   mean
  #   over time
  #   for some units
  plot_enveffect(val2 = val2, 
                 data = data, 
                 sam = sam, 
                 wei.fit = wei.fit, 
                 type = type, 
                 name = name, 
                 lump = lump,
                 npred = npred)
}
