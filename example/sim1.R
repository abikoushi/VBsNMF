library(NMF)
library(VBsNMF)
library(Matrix)
library(ggplot2)

set_data <- function(L, nrow, ncol, center=0, scale=1){
  Z <- matrix(rnorm(L*nrow,0,scale), nrow, L)
  Z <- sweep(Z,1,rowMeans(Z)-center)
  W <- matrix(rnorm(L*ncol,0,scale), L, ncol)
  W <- sweep(W,1,rowMeans(W)-center)
  Y <- matrix(rpois(nrow*ncol, exp(Z)%*%exp(W)), nrow, ncol)
  list(Y=Y,Z=Z,W=t(W))
}

#####
#VB vs. EM
set.seed(1); dat <- set_data(2, 100, 100)
system.time({
  out_vb <- VBNMF(dat$Y, rank=2, iter=250, prior_rate=1)
})
plot(out_vb$logprob[-1], type="l")

system.time({
  out_em <- em_nmf_pois(dat$Y, rank=2, iter=100, prior_rate=1)
})
plot(out_em$logprob, type="l")

sqrt(mean((dat$Y-basemean(out_vb)%*%coefmean(out_vb))^2))
sqrt(mean((dat$Y-out_em$Z%*%t(out_em$W))^2))

###
#missing
indna <- sample.int(length(dat$Y), size = 100)
Y <- dat$Y
Y[indna] <- NA

system.time({
  out <- vb_nmf_pois(dat$Y, rank=2, iter=1000, prior_rate=1)
})
plot(out$logprob)
Zhat <- basemean(out)
What <- coefmean(out)

fit <- Zhat%*%What

basecol <-rgb(0,0,0,0.1)
plot(dat$Y[-indna], fit[-indna], col=basecol, cex=0.5)
points(dat$Y[indna], fit[indna], col=rgb(1,0.5,0,0.3), pch=16)
abline(0,1, col="lightgrey", lty=2)
####

