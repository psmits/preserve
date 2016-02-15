# R code for JADE (Joint species-rank Abundance Distribution/Estimation) based on individual-based (abundance) and sampling-unit-based (incidence) data.

#' DetAbu(x, zero=FALSE) is a function of estimating detected species relative abundance.
#' @param x a vector of species abundance frequency
#' @param zero reserves zero frequency or not. Default is FALSE.
#' @return a numerical vector
DetAbu <- function(x, zero=FALSE){
  x <- unlist(x)
  n <- sum(x)  
  f1 <- sum(x==1)
  f2 <- sum(x==2)
  f3 <- sum(x==3)
  if(f2==0){
    f1 <- max(f1 - 1, 0)
    f2 <- 1
  }
  A1 <- f1 / n * ((n-1)*f1 / ((n-1)*f1 + 2*max(f2,1)))
  A2 <- f2 / choose(n, 2) * ((n-2)*f2 / ((n-2)*f2 + 3*max(f3,1)))^2
  if(zero==FALSE) x <- x[x>0]
  q.solve <- function(q){
    e <- A1 / sum(x/n*exp(-q*x))
    out <- sum((x/n * (1 - e * exp(-q*x)))^2) - sum(choose(x,2)/choose(n,2)) + A2
    abs(out)
  }
  q <- tryCatch(optimize(q.solve, c(0,1))$min, error = function(e) {1})
  e <- A1 / sum(x/n*exp(-q*x))
  o <- x/n * (1 - e * exp(-q*x))
  o
}

#' UndAbu(x) is a function of estimating undetected species relative abundance.
#' @param x a vector of species abundance frequency
#' @return a numerical vector
UndAbu <- function(x){
  x <- unlist(x)
  n <- sum(x)  
  f1 <- sum(x==1)
  f2 <- sum(x==2)
  f3 <- sum(x==3)
  f4 <- max(sum(x == 4), 1)
  f0.hat <- ceiling(ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2))  # Chao1 estimator
  if(f2==0){
    f1 <- max(f1 - 1, 0)
    f2 <- 1
  }
  A1 <- f1 / n * ((n-1)*f1 / ((n-1)*f1 + 2*max(f2,1)))
  A2 <- f2 / choose(n, 2) * ((n-2)*f2 / ((n-2)*f2 + 3*max(f3,1)))^2
  R <- A1^2/A2
  j <- 1:f0.hat
  f.solve <- function(x){ 
    out <- sum(x^j)^2 / sum((x^j)^2) - R
    abs(out)
  }
  b <-  tryCatch(optimize(f.solve, lower=(R-1)/(R+1), upper=1, tol=1e-5)$min, error = function(e) {(R-1)/(R+1)})
  a <- A1 / sum(b^j)
  p <- a * b^j
  if(f0.hat == 1) p <- A1
  p
}

#' DetInc(y, zero=FALSE) is a function of estimating detected species incidence probability.
#' @param y a vector of species incidence frequency
#' @param zero reserves zero frequency or not. Default is FALSE.
#' @return a numerical vector
DetInc <- function(y, zero=FALSE){
  y <- unlist(y)
  nT <- max(y)  
  y <- y[-1]
  Q1 <- sum(y==1)
  Q2 <- sum(y==2)
  Q3 <- sum(y==3)
  if(Q2==0){
    Q1 <- max(Q1 - 1, 0)
    Q2 <- 1
  }
  A1 <- Q1 / nT * ((nT-1)*Q1 / ((nT-1)*Q1 + 2*max(Q2,1)))
  A2 <- Q2 / choose(nT, 2) * ((nT-2)*Q2/((nT-2)*Q2 + 3*max(Q3,1)))^2
  if(zero==FALSE) y <- y[y>0]
  q.solve <- function(q){
    e <- A1 / sum(y/T*exp(-q*y))
    out <- sum((y/nT * (1 - e * exp(-q*y)))^2) - sum(choose(y,2)/choose(nT,2)) + A2
    abs(out)
  }
  q <- tryCatch(optimize(q.solve, c(0,1))$min, error = function(e) {1})
  e <- A1 / sum(y/nT*exp(-q*y))
  o <- y/nT * (1 - e * exp(-q*y))
  o
}

#' UndInc(y) is a function of estimating undetected species incidence probability.
#' @param y a vector of species incidence frequency.
#' @return a numerical vector
UndInc <- function(y){
  y <- unlist(y)
  nT <- max(y)  
  y <- y[-1]
  Q1 <- sum(y==1)
  Q2 <- sum(y==2)
  Q3 <- sum(y==3)
  Q4 <- max(sum(y == 4), 1)
  Q0.hat <- ceiling(ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2))  #Chao2 estimator
  if(Q2==0){
    Q1 <- max(Q1 - 1, 0)
    Q2 <- 1
  }
  A1 <- Q1 / nT * ((nT-1)*Q1 / ((nT-1)*Q1 + 2*max(Q2,1)))
  A2 <- Q2 / choose(nT, 2) * ((nT-2)*Q2/((nT-2)*Q2 + 3*max(Q3,1)))^2
  R <- A1^2/A2
  j <- 1:Q0.hat
  f.solve <- function(x){ 
    out <- sum(x^j)^2 / sum((x^j)^2) - R
    abs(out)
  }
  b <-  tryCatch(optimize(f.solve, lower=(R-1)/(R+1), upper=1, tol=1e-5)$min, error = function(e) {(R-1)/(R+1)})
  a <- A1 / sum(b^j)
  p <- a * b^j
  if(Q0.hat ==1) p <- A1
  p
}


#' SpecDist(x, datatype=¡¨abundance¡¨) is a function of estimating species rank distribution for abundance/incidence based data
#' @param x a vector of species abundance or incidence frequency. If datatype = "incidence", then the input format of first entry should be total number of sampling units, and followed by species incidence frequency.
#' @param datatype the data type (abundance or incidence) of input.
#' @return a data.frame object of RAD/RID
SpecDist <- function(x, datatype="abundance"){
  TYPE <- c("abundance", "incidence")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  
  if(datatype=="abundance"){  
    out <- rbind(data.frame("probability"=DetAbu(x, zero=TRUE), 
                            "method"="detected"),
                 data.frame("probability"=UndAbu(x), 
                            "method"="undetected"))
  }else if(datatype=="incidence"){
    out <- rbind(data.frame("probability"=DetInc(x, zero=TRUE), 
                            "method"="detected"),
                 data.frame("probability"=UndInc(x), 
                            "method"="undetected"))
  }
  out[order(-out[,1]),]
}

## Examples
#x <- read.csv("https://raw.githubusercontent.com/JohnsonHsieh/JADE/master/sandbox/ex-abun.csv")
#out1 <- SpecDist(x, "abundance")
#out1 <- subset(out1, probability>0)
#plot(out1$probability, col=ifelse(out1$method=="detected",1,2), xlab="Species Rank", ylab="Relative Abundance")
# 
#
#y <- read.csv("http://raw.githubusercontent.com/JohnsonHsieh/JADE/master/sandbox/ex-inci.csv")
#out2 <- SpecDist(y, "incidence")
#out2 <- subset(out2, probability>0)
#plot(out2$probability, col=ifelse(out2$method=="detected",1,2), xlab="Species Rank", ylab="Incidence Probability")
# 
