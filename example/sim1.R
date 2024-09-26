library(NMF)
library(VBsNMF)
library(rbenchmark)
library(bench)

set_data <- function(L, nrow, ncol, scale){
  Z <- matrix(rgamma(L*nrow,1,1), nrow, L)
  #Z <- Z/rowSums(Z)
  W <- matrix(rgamma(L*ncol,1,1), L, ncol)
  #W <- W/rowSums(W)
  Y <- matrix(rpois(nrow*ncol, (Z%*%W)*scale), nrow, ncol)
  list(Y=Y,Z=Z,W=t(W))
}

set.seed(1010); dat <- set_data(2, 100, 100, scale = 1)

mean(dat$Y==0)
Y <- as(dat$Y, "TsparseMatrix")

system.time({
  out <- vb_nmf_pois(Y, rank=2, iter=10, prior_rate=0.1)
})

system.time({
  out2 <- nmf(dat$Y, rank=2, maxIter=10, eps=0)
})

system.time({
  out_s <- svb_nmf_pois(Y, rank=2, b_size = 100)
})


plot(out_s$logprob, type = "l", lty=1)

Zhat <- basemean(out)
What <- coefmean(out)

plot(dat$Y, basis(out2)%*%coef(out2), col=rgb(1,0,0,0.1), pch=2)
points(dat$Y, Zhat%*%What, col=rgb(0,0,1,0.1))

sqrt(mean((dat$Y-basis(out2)%*%coef(out2))^2))
sqrt(mean((dat$Y-Zhat%*%What)^2))

Zhat_s <- basemean(out_s)
What_s <- coefmean(out_s)

plot(Zhat_s%*%What_s, Y, col=rgb(0,0,0,0.1))
points(dat$Y, Zhat%*%What, col=rgb(0,0,1,0.1),pch=5)
points(dat$Y, basis(out2)%*%coef(out2), col=rgb(1,0,0,0.1), pch=2)

sqrt(mean((dat$Y-basis(out2)%*%coef(out2))^2))
sqrt(mean((dat$Y-Zhat%*%What)^2))
sqrt(mean((dat$Y-Zhat_s%*%What_s)^2))

