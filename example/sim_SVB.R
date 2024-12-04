library(Matrix)
library(ggplot2)
library(VBsNMF)

set_data <- function(L, nrow, ncol, center=0, scale=1){
  Z <- matrix(rnorm(L*nrow,0,scale), nrow, L)
  Z <- sweep(Z,1,rowMeans(Z)-center)
  W <- matrix(rnorm(L*ncol,0,scale), L, ncol)
  W <- sweep(W,1,rowMeans(W)-center)
  Y <- matrix(rpois(nrow*ncol, exp(Z)%*%exp(W)), nrow, ncol)
  list(Y=Y,Z=Z,W=t(W))
}

set.seed(1234); dat <- set_data(2, 101, 105)
length(dat$Y)
Y <- as(dat$Y, "TsparseMatrix")
system.time({
  out <- VBsNMF:::VBNMF(Y, rank = 2)
})

plot(out$logprob, type="l")

Zhat <- basemean(out)
What <- coefmean(out)

plot(dat$Y, Zhat%*%What, pch=".")
abline(0,1)
####

writeMM(Y,"test.mtx")
hist(dat$Y, breaks = "FD")

system.time({
out_s <- VBsNMF:::SVBNMF("test.mtx",
                         rank = 2,
                         b_size = 100,
                         subiter = 5,
                         forgetting = 0.9, delay = 1,
                         n_epochs = 200)
})

plot(out_s$logprob[-1], type="l")

colSums(basemean(out))
rowSums(coefmean(out))

Zhat <- basemean(out_s)
What <- coefmean(out_s)
fit <- Zhat%*%What

colSums(Zhat)
rowSums(What)

ggplot(data=NULL, aes(x=c(dat$Y), y=c(fit)))+
  geom_abline(intercept = 0, slope = 1, col="grey", linetype=1) +
  geom_bin2d(aes(fill = after_stat(log2(count))))+
  theme_classic(16)+labs(x="observerd", y="fitted")

# ggplot(data=NULL, aes(x=c(dat$Y), y=c(fit)))+
#   geom_abline(intercept = 0, slope = 1, col="royalblue", linetype=1) +
#   geom_point(size=0.1, alpha=0.5)+
#   theme_classic(16)+labs(x="observerd", y="fitted")
