library(NMF)
library(VBsNMF)

set_data <- function(L, nrow, ncol, center=0, scale=1){
  Z <- matrix(rnorm(L*nrow,0,scale), nrow, L)
  Z <- sweep(Z,1,rowMeans(Z)-center)
  W <- matrix(rnorm(L*ncol,0,scale), L, ncol)
  W <- sweep(W,1,rowMeans(W)-center)
  Y <- matrix(rpois(nrow*ncol, exp(Z)%*%exp(W)), nrow, ncol)
  list(Y=Y,Z=Z,W=t(W))
}

set.seed(1111); dat <- set_data(2, 200, 200)


indna <- sample.int(length(dat$Y), size = 10)
Y <- dat$Y
Y[indna] <- NA
#mean(dat$Y==0)
Y <- as(Y, "TsparseMatrix")
hist(dat$Y, breaks = "FD")


system.time({
  out <- vb_nmf_pois(Y, rank=2, iter=100, prior_rate=1)
})


plot(out$logprob)

Zhat <- basemean(out)
What <- coefmean(out)

fit <- (Zhat%*%What)

basecol <-rgb(0,0,0,0.1)
plot(dat$Y[-indna], fit[-indna], col=basecol, cex=0.5)
points(dat$Y[indna], fit[indna], col=rgb(1,0.5,0,0.5), pch=16)
abline(0,1, col="royalblue")

system.time({
  out <- vb_nmf_pois(Y, rank=2, iter=10, prior_rate=0.1)
})

system.time({
  out2 <- nmf(dat$Y, rank=2, maxIter=10, eps=0)
})

system.time({
out_s <- svb_nmf_pois(Y, rank=2, b_size = 100, forgetting = 0.9, delay = 1.5, n_epochs = 100)
})

plot(out_s$logprob)
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

