# basic posterior predictive checks based on bayesplot
plot_ppcbasics <- function(wei.fit, data, lump, mvgr, cn, name) {
  duration <- c(data$dur)

  pm <- ppc_stat(duration, wei.fit$y_tilde, 'mean')
  pm <- pm + labs(x = 'Duration (geological stages)',
                  y = 'Frequency')
  ggsave(pm, filename = paste0('../doc/figure/ppc_mean_',
                               name, '.pdf'),
         width = 6, height = 5, dpi = 600)

  pmg <- ppc_stat_grouped(duration, wei.fit$y_tilde, 
                          group = mvgr, 'mean',
                          facet_args = list(strip.position = 'bottom'))
  pmg <- pmg + labs(x = 'Duration (geological stages)',
                    y = 'Frequency')
  pmg <- pmg + theme(strip.text = element_text(size = 10))
  ggsave(pmg, filename = paste0('../doc/figure/ppc_mean_group_',
                                name, '.pdf'),
         width = 10, height = 8, dpi = 600)

  pe <- ppc_stat(duration, wei.fit$y_tilde, 'median')
  pe <- pe + labs(x = 'Duration (geological stages)',
                  y = 'Frequency')
  ggsave(pe, filename = paste0('../doc/figure/ppc_median_',
                               name, '.pdf'),
         width = 6, height = 5, dpi = 600)


  peg <- ppc_stat_grouped(duration, wei.fit$y_tilde, 
                          group = mvgr, 'median',
                          facet_args = list(strip.position = 'bottom'))
  peg <- peg + labs(x = 'Duration (geological stages)',
                    y = 'Frequency')
  peg <- peg + theme(strip.text = element_text(size = 10))
  ggsave(peg, filename = paste0('../doc/figure/ppc_med_group_',
                                name, '.pdf'),
         width = 10, height = 8, dpi = 600)

  pd <- ppc_dens_overlay(duration, wei.fit$y_tilde[1:100, ])
  pd <- pd + labs(x = 'Duration (geological stages)',
                  y = 'Density')
  ggsave(pd, filename = paste0('../doc/figure/ppc_dens_',
                               name, '.pdf'),
         width = 6, height = 5, dpi = 600)
  pd <- pd + coord_cartesian(xlim = c(0, 30))
  ggsave(pd, filename = paste0('../doc/figure/ppc_dens_zoom_',
                               name, '.pdf'),
         width = 6, height = 5, dpi = 600)

  ps <- ppc_error_scatter_avg(duration, wei.fit$y_tilde)
  ggsave(ps, filename = paste0('../doc/figure/ppc_err_scatter_avg_',
                               name, '.pdf'),
         width = 6, height = 5, dpi = 600)

  pe <- ppc_ecdf_overlay(duration, wei.fit$y_tilde[1:100, ])
  ggsave(pe, filename = paste0('../doc/figure/ppc_ecdf_',
                               name, '.pdf'),
         width = 6, height = 5, dpi = 600)
  pe <- pe + coord_cartesian(xlim = c(0, 30))
  ggsave(pe, filename = paste0('../doc/figure/ppc_ecdf_zoom_',
                               name, '.pdf'),
         width = 6, height = 5, dpi = 600)
}




# these are possible values/partially counter factual statements
# environmental effect
# get the "mustache" for each data point
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

# get the mean mustache
quad.mean <- function(x, y, mcoef, npred, sam) {
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

# plots of the effects of environmental preference 
#   mean
#   over time
#   for some units
# this is called for side-effects
plot_enveffect <- function(val2, # quantile
                           data, # data fit by model
                           sam, # which samples?
                           wei.fit, # model posterior
                           type, # type of analysis
                           name, # name of model
                           lump, # translating times
                           npred) # number of predictors
{
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

    mcoef <- colMeans(wei.fit$mu_prior[sam, ])[1:5]
    meanquad <- data.frame(env = val, 
                           resp = quad.mean(val, y = val2[zz], 
                                            mcoef, npred = npred,
                                            sam = sam))

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




# plots of the covariate effects estimates over time
plot_coveffect <- function(wei.fit, npred, data, lump, name) {
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
                                          scales = 'free_y', 
                                          switch = 'y',
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
}




# extract covariance/correlation matrices for inspection
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
