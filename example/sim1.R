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

set.seed(1234); dat <- set_data(2, 99, 205)
#mean(dat$Y==0)
Y <- as(dat$Y, "TsparseMatrix")
writeMM(Y,"test.mtx")
out_s <- VBsNMF:::svb_nmf_pois_mtx("test.mtx", 2, b_size = 1000,
                                   forgetting = 0.9, delay = 1,
                                   n_epochs = 100)
hist(dat$Y, breaks = "FD")
plot(out_s$logprob, type="l")
Zhat <- basemean(out_s)
What <- coefmean(out_s)

fit <- Zhat%*%What
length(dat$Y)
length(fit)
basecol <-rgb(0,0,0,0.1)
ggplot(data=NULL, aes(x=c(dat$Y), y=c(fit)))+
  geom_abline(intercept = 0, slope = 1, col="lightgrey", linetype=2) +
  geom_bin2d(aes(fill=after_stat(log10(count))))+
  theme_classic(16)+labs(x="observerd", y="fitted")

#####
system.time({
  out <- vb_nmf_pois(dat$Y, rank=2, iter=100, prior_rate=1)
})

system.time({
  model <- RcppML::nmf(dat$Y, k = 2, maxit = 100, tol=0)
})
###
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

# set.seed(1234);
# dat <- rsparsematrix(1000,1000,0.9,
#                      rand.x = function(n){rpois(n,10)+1})

set.seed(1); dat <- set_data(2, 100, 100)
system.time({
  out <- VBsNMF::vb_nmf_pois(dat$Y, rank=2,
                             prior_shape=1, prior_rate=1, iter=1000)
})

plot(out$logprob)
system.time({
  out2 <- NMF::nmf(dat$Y, rank=2, maxIter=1000, eps=-1)
})

Zhat <- basemean(out)
What <- coefmean(out)

# Zhat_s <- basemean(out_s)
# What_s <- coefmean(out_s)
#plot(Y, Zhat_s%*%What_s,  col=rgb(0,0,0,0.1))
plot(dat$Y, Zhat%*%What, col=rgb(0,0,1,0.1), pch=5)
points(dat$Y, basis(out2)%*%coef(out2), col=rgb(1,0,0,0.1), pch=2)
abline(0,1,lty=2)

sqrt(mean((dat$Y-basis(out2)%*%coef(out2))^2))
sqrt(mean((dat$Y-Zhat%*%What)^2))
#sqrt(mean((dat$Y-Zhat_s%*%What_s)^2))

