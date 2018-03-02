# start the whole plotting sequence
#########################
source('../R/multiplot.r')
source('../R/borrow_plotcorr.r')
posterior.plots <- function(data, wei.fit, npred, name = 'cweib') {
  # data setup
  # HERE
  coh <- c(data$cohort)
  rage <- c(data$occupy)
  envs <- c(data$env)
  size <- c(data$size)
  duration <- c(data$dur)


  # data distributions would be cool, specifically in the case of the imputation

  ss <- sample(length(wei.fit$lp__), 12)
  impres <- data.frame(melt(wei.fit$samp[ss, ]), rep(data$inclusion, each = 12))
  names(impres) <- c('sim', 'var', 'value', 'inc')
  imp.gg <- ggplot(impres, aes(x = value, fill = factor(inc)))
  imp.gg <- imp.gg + geom_histogram(bins = 10)
  imp.gg <- imp.gg + facet_wrap( ~ sim) 
  imp.gg <- imp.gg + scale_fill_manual(values = c('grey', 'black'),
                                       name = 'Data source',
                                       labels = c('Imputed', 'Observed'))
  imp.gg <- imp.gg + labs(x = 'Gap statistic', y = 'Count')
  imp.gg <- imp.gg + theme(strip.text = element_blank(),
                           strip.background = element_blank(),
                           legend.title = element_text(size = 10),
                           legend.text = element_text(size = 8))
  ggsave(imp.gg, filename = paste0('../doc/figure/imputation_compare_',
                                   name, '.pdf'),
         width = 8, height = 4, dpi = 600)



  # lets make survival curves
  # non-parametric estimate of the survival curve just given durations
  #   some durations are left-censored
  condition <- (data$censored == 0) * 1
  emp.surv <- survfit(Surv(time = duration, event = condition, type = 'left') ~ 1)
  emp.surv <- data.frame(time = emp.surv$time, surv = emp.surv$surv, 
                         lower = emp.surv$lower, upper = emp.surv$upper)

  # i have two ways of estimating the survival functions from posterior
  #   actual S(t) estimates
  #   estimated durations under discrete Weibull
  wr <- wei.fit$y_tilde[sample(nrow(wei.fit$y_tilde), 100), ]
  wei.surv <- apply(wr, 1, function(x) survfit(Surv(x) ~ 1))
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

  # in b&w
  surv.plot <- ggplot(emp.surv, aes(x = time, y = surv))
  surv.plot <- surv.plot + geom_line(data = sim.surv, 
                                     aes(x = time, y = surv, group = group),
                                     colour = 'black', alpha = 0.05)
  surv.plot <- surv.plot + geom_line(colour = 'black', size = 1.2)
  surv.plot <- surv.plot + coord_cartesian(xlim = c(-0.5, max(duration) + 2))
  surv.plot <- surv.plot + labs(x = 'Duration (t)', 
                                y = 'Pr(t < T)')
  surv.plot <- surv.plot + theme(axis.title = element_text(size = 25))
  surv.plot <- surv.plot + scale_y_continuous(trans=log10_trans(),
                                              breaks = c(0.01, 0.1, 0.5, 1))
  ggsave(surv.plot, filename = paste0('../doc/figure/survival_curves_bw_',
                                      name, '.pdf'),
         width = 6, height = 5, dpi = 600)


  # posterior predictive point checks
  quant <- laply(wr, function(x) quantile(x, seq(0.1, 0.9, by = 0.05)))
  qudur <- quantile(duration, seq(0.1, 0.9, by = 0.05))
  qp <- colSums(t(apply(quant, 1, function(x) x > qudur))) / nrow(quant)

  point.check <- data.frame(quan = seq(0.1, 0.9, by = 0.05), percent = qp)
  point.plot <- ggplot(point.check, aes(x = quan, y = percent))
  point.plot <- point.plot + geom_point()
  point.plot <- point.plot + labs(x = 'quantile of observed duration',
                                  y = 'P(quantile of simulated > 
                                  quantile of observed)')
  ggsave(point.plot, filename = paste0('../doc/figure/quantile_', name, '.pdf'),
         width = 6, height = 5, dpi = 600)


  est.shotgun <- data.frame(obs = duration, 
                            sim = colMeans(wei.fit$hold))

  shot.plot <- ggplot(est.shotgun, aes(x = sim, y = obs))
  shot.plot <- shot.plot + stat_function(fun = function(x) x, 
                                         lty = 'dashed', 
                                         colour = 'darkgrey')
  shot.plot <- shot.plot + geom_point(alpha = 0.5)
  shot.plot <- shot.plot + labs(x = 'Estimated duration approx. (t)',
                                y = 'Observed duration (t)')
  #shot.plot <- shot.plot + labs(x = expression(tilde(sigma)),
  #                              y = 'Observed duration (t)')
  ggsave(shot.plot, filename = paste0('../doc/figure/shotgun_', name, '.pdf'),
         width = 6, height = 5, dpi = 600)

  # quality of fit is medium, though a lot is captured


  # make plot of correlation and covariance matrices
  # row is sample
  # dim 2 is row
  # dim 3 is col
  get.covcor <- function(stanfit, npred) {
    cor.median <- matrix(, ncol = npred, nrow = npred)
    cor.mean <- matrix(, ncol = npred, nrow = npred)
    cor.10 <- matrix(, ncol = npred, nrow = npred)
    cor.90 <- matrix(, ncol = npred, nrow = npred)
    cov.median <- matrix(, ncol = npred, nrow = npred)
    cov.mean <- matrix(, ncol = npred, nrow = npred)
    cov.10 <- matrix(, ncol = npred, nrow = npred)
    cov.90 <- matrix(, ncol = npred, nrow = npred)
    for(ii in seq(npred)) {
      for(jj in seq(npred)) {
        cor.median[jj, ii] <- median(stanfit$Omega[, jj, ii])
        cor.mean[jj, ii] <- mean(stanfit$Omega[, jj, ii])
        cor.10[jj, ii] <- quantile(stanfit$Omega[, jj, ii], probs = .1)
        cor.90[jj, ii] <- quantile(stanfit$Omega[, jj, ii], probs = .9)
        cov.median[jj, ii] <- median(stanfit$Sigma[, jj, ii])
        cov.mean[jj, ii] <- mean(stanfit$Sigma[, jj, ii])
        cov.10[jj, ii] <- quantile(stanfit$Sigma[, jj, ii], probs = .1)
        cov.90[jj, ii] <- quantile(stanfit$Sigma[, jj, ii], probs = .9)
      }
    }

    out <- list(cor.median, cor.mean, cor.10, cor.90, 
                cov.median, cov.mean, cov.10, cov.90)
    if(npred == 6) {
      out <- llply(out, function(x) {
                     nn <- c('i', 'r', 'e', 'e2', 'rxe', 'm', 's')
                     nn <- nn[seq(npred)]
                     rownames(x) <- colnames(x) <- nn
                     x})
    } else if(npred == 5) {
      out <- llply(out, function(x) {
                     nn <- c('i', 'r', 'e', 'e2', 'm', 's')
                     nn <- nn[seq(npred)]
                     rownames(x) <- colnames(x) <- nn
                     x})
    }
    out
  }
  wei.covcor <- get.covcor(wei.fit, npred)

  # just for the weibull
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

  # effect of all the covariates
  # mean cohort effects
  efmu <- colMeans(wei.fit$mu_prior)
  efmurange <- apply(wei.fit$mu_prior, 2, function(x) 
                     quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9)))

  # effects for each cohort
  efbeta <- colMeans(wei.fit$beta)
  efbetarange <- efbetaprob <- list()
  for(ii in seq(data$O)) {
    efbetarange[[ii]] <- apply(wei.fit$beta[, ii, ], 2, function(x) 
                               quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9)))
    efbetaprob[[ii]] <- apply(wei.fit$beta[, ii, ], 2, function(x) 
                              sum(x > 0) / length(x))
  }

  ef.df <- t(rbind(mean = efmu, efmurange, pred = seq(npred)))
  ef.df <- cbind(rbind(ef.df, ef.df), 
                 time = c(rep(1, times = npred), rep(data$O, times = npred)))
  ef.df <- data.frame(ef.df)
  #ef.df$time

  efbeta.h <- Map(function(x) t(rbind(mean = efbeta[x, ], efbetarange[[x]])), 
                  seq(data$O))
  efbeta.h <- Map(function(x) data.frame(efbeta.h[[x]], 
                                         time = x, 
                                         pred = seq(npred)), seq(data$O))
  efbeta.df <- Reduce(rbind, efbeta.h)
  if(npred == 6) {
    too <- c('intensity', 'range', 'env_pref', 'env_curv', 
             'rxe', 'size', 'delta')[seq(npred)]
  } else if(npred == 5) {
    too <- c('intensity', 'range', 'env_pref', 'env_curv', 
             'size', 'delta')[seq(npred)]
  }
  #too <- c('beta^0', 'beta^r', 'beta^v', 'beta^v^2', 'beta^m', 
  #         'delta')[seq(npred)]
  efbeta.df$pred <- mapvalues(efbeta.df$pred, 
                              from = seq(npred), 
                              to = too)
  ef.df$pred <- mapvalues(ef.df$pred, 
                          from = seq(npred), 
                          to = too)

  efbeta.df$pred <- factor(efbeta.df$pred, levels = too)
  ef.df$pred <- factor(ef.df$pred, levels = too)

  efbeta.df$time <- mapvalues(efbeta.df$time, seq(33), lump[5:(5+33-1), 3])
  ef.df$time <- mapvalues(ef.df$time, seq(33), lump[5:(5+33-1), 3])


  efbeta.plot <- ggplot(efbeta.df, aes(x = time, y = X50.))
  efbeta.plot <- efbeta.plot + geom_pointrange(mapping = aes(ymin = X10., 
                                                             ymax = X90.),
                                               fatten = 2)
  efbeta.plot <- efbeta.plot + facet_grid(pred ~ .,
                                          scales = 'free_y', switch = 'y',
                                          labeller = label_parsed)
  efbeta.plot <- efbeta.plot + geom_ribbon(data = ef.df, 
                                           mapping = aes(ymin = X10.,
                                                         ymax = X90.),
                                           alpha = 0.2)
  efbeta.plot <- efbeta.plot + geom_line(data = ef.df, 
                                         mapping = aes(y = X50.),
                                         alpha = 0.5)
  efbeta.plot <- efbeta.plot + labs(x = 'Time (My)', y = 'Effect estimate for...')
  efbeta.plot <- efbeta.plot + scale_x_reverse()
  ggsave(efbeta.plot, filename = paste0('../doc/figure/cohort_series_', 
                                        name, '.pdf'),
         width = 7.5, height = 10, dpi = 600)
  ggsave(efbeta.plot, filename = paste0('../doc/figure/cohort_series_wide_', 
                                        name, '.pdf'),
         width = 12, height = 10, dpi = 600)



  # these are possible values/partially counter factual statements
  # environmental effect
  quad <- function(x, y, sam, npred) {
    bet <- wei.fit$mu_prior[sam, 1]
    if(npred == 6) {
      bet <- -(bet + (wei.fit$mu_prior[sam, 2] * y) +
               (wei.fit$mu_prior[sam, 3] * x) + 
               (wei.fit$mu_prior[sam, 4] * x^2) +
               (wei.fit$mu_prior[sam, 5] * (x * y)))
    } else if (npred == 5) {
      bet <- -(bet + (wei.fit$mu_prior[sam, 2] * y) +
               (wei.fit$mu_prior[sam, 3] * x) + 
               (wei.fit$mu_prior[sam, 4] * x^2))
    }
    if(is.null(wei.fit$alpha)) {
      o <- bet
    } else {
      o <- bet / wei.fit$alpha[sam]  # depends on if alpha varies by cohort
    }
    o
  }
  quad.mean <- function(x, y, mcoef, npred) {
    #if(!(best %in% 1:2)) {
    #  -(mcoef[1] + (mcoef[2] * x) + (mcoef[3] * x^2)) / exp(mean(wei.fit$alpha_mu[sam]))
    #} else {
    if(npred == 6) {
      bet <- -(mcoef[1] + (mcoef[2] * y) + (mcoef[3] * x) + (mcoef[4] * x^2) + 
               (mcoef[5] * (x * y)))
    } else if(npred == 5) {
      bet <- -(mcoef[1] + (mcoef[2] * y) + (mcoef[3] * x) + (mcoef[4] * x^2))
    }
    if(is.null(wei.fit$alpha)) {
      o <- bet
    } else {
      o <- bet / mean(wei.fit$alpha[sam])
    }
    o
  }

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
  for(zz in seq(length(val2))) {
    #sam <- sample(nrow(wei.fit$lp__), 1000)

    # HERE
    env.d <- c(data$env)
    val <- seq(from = min(env.d), to = max(env.d), by = 0.01)
    quadval <- list()
    for(ii in seq(length(sam))) {
      quadval[[ii]] <- data.frame(env = val, 
                                  resp = quad(val, y = val2[zz], 
                                              sam[ii], npred = npred), 
                                  sim = ii)
    }
    quadframe <- Reduce(rbind, quadval)

    mcoef <- colMeans(wei.fit$mu_prior[sam, ])[c(1, 2, 3, 4, 5, 6)]
    meanquad <- data.frame(env = val, 
                           resp = quad.mean(val, y = val2[zz], 
                                            mcoef, npred = npred))

    # add rug showing observed
    #   this addition would overpower the big, by cohort graph
    # HERE
    env.obs <- data.frame(env = data$env)

    mustache <- ggplot(quadframe, aes(x = env, y = resp, group = sim))
    mustache <- mustache + geom_line(alpha = 1 / 100)
    mustache <- mustache + geom_line(data = meanquad,
                                     mapping = aes(group = NULL),
                                     colour = 'blue')
    mustache <- mustache + geom_rug(data = env.obs,
                                    mapping = aes(x = env, y = NULL, group = NULL),
                                    sides = 'b', alpha = 0.05)
    mustache <- mustache + labs(x = 'Environmental preference\n(open-ocean <--> epicontinental)', 
                                y = 'log(approx. expected duration in t)')
    #y = expression(paste('log(', sigma, ')')))
    mustache <- mustache + theme(axis.title.x = element_text(hjust = 0.5))
    ggsave(mustache, filename = paste0('../doc/figure/env_effect_', 
                                       type[zz], '_', name, '.pdf'),
           width = 6, height = 5, dpi = 600)



    # by cohort
    sam2 <- sam[1:100]
    if(npred == 6) {
      bet.coh <- wei.fit$beta[sam2, , c(1, 2, 3, 4, 5)]
    } else if(npred == 5) {
      bet.coh <- wei.fit$beta[sam2, , c(1, 2, 3, 4)]
    }
    #if(!(best %in% 1:2)) {
    #  alp.coh <- apply(wei.fit$alpha_cohort[sam, ], 2, function(x) 
    #                   x + wei.fit$alpha_mu[sam])
    #} else {
    alp.coh <- wei.fit$alpha_trans[sam2]
    #}
    val <- seq(from = min(env.d), to = max(env.d), by = 0.01)
    if(npred == 6) {
      dat <- cbind(1, val2[zz], val, val^2, val * val2[zz])
    } else if(npred == 5) {
      dat <- cbind(1, val2[zz], val, val^2)
    }

    coh.est <- list()
    for(ii in seq(data$O)) {
      h <- list()
      for(jj in seq(length(val))) {
        if(is.null(alp.coh)) {
          h[[jj]] <- -(bet.coh[, ii, ] %*% dat[jj, ])
        } else {
          h[[jj]] <- -(bet.coh[, ii, ] %*% dat[jj, ]) / exp(alp.coh[ii])
        }
      }
      coh.est[[ii]] <- h
    }

    # massage into shape
    #   val, resp (V2), sim, coh
    stg.name <- as.character(lump[5:(5+33-1), 4])
    stg.name <- Reduce(c, Map(function(x, y) paste0(y, '. ', x), 
                              stg.name, seq(length(stg.name))))
    coh.map <- list()
    for(jj in seq(data$O)) {
      h <- Map(function(x, y) {
                 cbind(val = x, resp = coh.est[[jj]][[y]], sim = seq(100))}, 
                 x = val, y = seq(length(val)))
      h <- Reduce(rbind, h)
      coh.map[[jj]] <- data.frame(h, coh = stg.name[jj])
    }
    coh.df <- Reduce(rbind, coh.map)

    coh.df.short <- coh.df[coh.df$coh %in% c('14. Emsian', '15. Eifelian', 
                                             '16. Givetian', '17. Frasnian'), ]

    cohmust <- ggplot(coh.df, aes(x = val, y = V2, group = sim))
    cohmust <- cohmust + geom_line(data = meanquad,
                                   mapping = aes(x = env,
                                                 y = resp,
                                                 group = NULL),
                                   colour = 'black', size = 1.5)
    cohmust <- cohmust + geom_line(alpha = 1 / 10, colour = 'blue')
    cohmust <- cohmust + facet_wrap(~ coh, strip.position = 'bottom', ncol = 7)
    cohmust <- cohmust + theme(axis.text = element_text(size = 8),
                               strip.text = element_text(size = 8))
    cohmust <- cohmust + labs(x = 'Environmental preference (v)',
                              y = 'log(approx. expected duration in t)')
    #y = expression(paste('log(', sigma, ')')))
    ggsave(cohmust, 
           filename = paste0('../doc/figure/env_cohort_', type[zz], 
                             '_', name, '.pdf'),
           width = 7.5, height = 8, dpi = 600)
    ggsave(cohmust, 
           filename = paste0('../doc/figure/env_cohort_wide_', type[zz], 
                             '_', name, '.pdf'),
           width = 9.5, height = 8, dpi = 600)
    cohmust.short <- cohmust %+% coh.df.short
    cohmust.short <- cohmust.short + theme(strip.text = element_text(size = 12))
    ggsave(cohmust.short, 
           filename = paste0('../doc/figure/env_cohort_short_', type[zz], 
                             '_', name, '.pdf'),
           width = 10, height = 5, dpi = 600)
  }
}
