library(plyr)

load('../data/epi_over_off.rdata')

cor.mean <- matrix(, ncol = 5, nrow = 5)
for(ii in seq(5)) {
  for(jj in seq(5)) {
    cor.mean[jj, ii] <- mean(wei.fit$Omega[, jj, ii])
  }
}
rownames(cor.mean) <- c('i', 'r', 'e', 'e2', 'm')
colnames(cor.mean) <- c('i', 'r', 'e', 'e2', 'm')

ll <- dim(wei.fit$Omega)[1]
p.cor.ir <- sum(wei.fit$Omega[, 1, 2] < 0) / ll
mean(wei.fit$Omega[, 1, 2])
sd(wei.fit$Omega[, 1, 2])
p.cor.ie2 <- sum(wei.fit$Omega[, 1, 4] < 0) / ll
mean(wei.fit$Omega[, 1, 4])
sd(wei.fit$Omega[, 1, 4])
p.cor.re2 <- sum(wei.fit$Omega[, 2, 4] > 0) / ll
mean(wei.fit$Omega[, 2, 4])
sd(wei.fit$Omega[, 2, 4])

mean(wei.fit$mu_prior[, 5])
sd(wei.fit$mu_prior[, 5])
sum(wei.fit$mu_prior[, 5] < 0) / ll

mean(wei.fit$mu_prior[, 3])
sd(wei.fit$mu_prior[, 3])
sum(wei.fit$mu_prior[, 3] < 0) / ll
