library(NMF)

set_data <- function(L, nrow, ncol, scale){
  Z <- matrix(rgamma(L*nrow,1,1),nrow,L)
  #Z <- Z/rowSums(Z)
  W <- matrix(rgamma(L*ncol,1,1),L,ncol)
  #W <- W/rowSums(W)
  Y <- matrix(rpois(nrow*ncol, (Z%*%W)*scale), nrow, ncol)
  list(Y=Y,Z=Z,W=t(W))
}

set.seed(999); dat <- set_data(2, 10000, 10000, scale = 1)

hist(dat$Y, breaks = "FD")
system.time({
  out <- vb_nmf_pois(dat$Y, L=2, iter=10)
})

system.time({
  out2 <- nmf(dat$Y, rank=2, maxIter=10, eps=0)
})

Zhat <- with(out, sweep(shape_row,1,c(rate_row),"/"))
What <- with(out, sweep(shape_col,1,c(rate_col),"/"))

plot(dat$Y, basis(out2)%*%coef(out2), col=rgb(1,0,0,0.2),pch=2,xlim=c(0,300),ylim=c(0,300))
points(dat$Y, Zhat%*%t(What), col=rgb(0,0,1,0.2),xlim=c(0,300),ylim=c(0,300))

head(out$logprob)
plot(rowSums(out$logprob)[-1], type = "l", lty=1)

