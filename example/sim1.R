library(NMF)
library(VBsNMF)
library(rbenchmark)
library(bench)

curve(learning_rate(x,1,0.9),0,100)
curve(learning_rate(x,1,0.99),add=TRUE)
curve(learning_rate(x,1,0.55),add=TRUE)

set_data <- function(L, nrow, ncol, center=0, scale=1){
  Z <- matrix(rnorm(L*nrow,0,scale), nrow, L)
  Z <- sweep(Z,1,rowMeans(Z)-center)
  W <- matrix(rnorm(L*ncol,0,scale), L, ncol)
  W <- sweep(W,1,rowMeans(W)-center)
  Y <- matrix(rpois(nrow*ncol, exp(Z)%*%exp(W)), nrow, ncol)
  list(Y=Y,Z=Z,W=t(W))
}

set.seed(1111); dat <- set_data(2, 100, 100)

#mean(dat$Y==0)
Y <- as(dat$Y, "TsparseMatrix")
hist(dat$Y, breaks = "FD")

N1 <- length(Y@x)

system.time({
  out <- vb_nmf_pois(Y, rank=2, iter=10, prior_rate=0.1)
})

system.time({
  out2 <- nmf(dat$Y, rank=2, maxIter=10, eps=0)
})

system.time({
out_s <- svb_nmf_pois(Y, rank=2, b_size = 10, forgetting = 0.8, delay = 2)
})

Zhat <- basemean(out)
What <- coefmean(out)

Zhat_s <- basemean(out_s)
What_s <- coefmean(out_s)
plot(Y, Zhat_s%*%What_s,  col=rgb(0,0,0,0.1))
points(dat$Y, Zhat%*%What, col=rgb(0,0,1,0.1), pch=5)
points(dat$Y, basis(out2)%*%coef(out2), col=rgb(1,0,0,0.1), pch=2)
abline(0,1,lty=2)

sqrt(mean((dat$Y-basis(out2)%*%coef(out2))^2))
sqrt(mean((dat$Y-Zhat%*%What)^2))
sqrt(mean((dat$Y-Zhat_s%*%What_s)^2))

