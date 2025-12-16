library(Matrix)
library(ggplot2)
library(VBsNMF)
library(rbenchmark)
library(pryr)
library(bench)

set_data <- function(L, nrow, ncol, center=0, scale=1){
  Z <- matrix(rnorm(L*nrow,0,scale), nrow, L)
  Z <- sweep(Z,1,rowMeans(Z)-center)
  W <- matrix(rnorm(L*ncol,0,scale), L, ncol)
  W <- sweep(W,1,rowMeans(W)-center)
  Y <- matrix(rpois(nrow*ncol, exp(Z)%*%exp(W)), nrow, ncol)
  list(Y=Y, Z=Z, W=t(W))
}

set.seed(1234); dat <- set_data(2, 101, 201)
length(dat$Y)
Y <- as(dat$Y, "TsparseMatrix")
system.time({
  out <- VBsNMF:::VBNMF(Y, rank = 2)
})

plot(out$logprob, type="l")

fit <-  basemean(out)%*%coefmean(out)


####

writeMM(Y,"test.mtx")
hist(dat$Y, breaks = "FD")
length(dat$Y)
bench::mark({
out_s <- VBsNMF:::SVBNMF_mtx("test.mtx",
                         rank = 2,
                         lr_param = c(15,0.8),
                         b_size = 1000,
                         subiter = 1,
                         n_epochs = 200)
}, iterations = 1)

plot(out_s$logprob[-1], type="l")

Zhat <- basemean(out_s)
What <- coefmean(out_s)
fit_s <- Zhat%*%What

colSums(Zhat)
rowSums(What)

ggplot(data=NULL, aes(x=c(dat$Y), y=c(fit_s)))+
  geom_abline(intercept = 0, slope = 1, alpha=0.2, linetype=1) +
  geom_bin2d(aes(fill = after_stat(log2(count))), binwidth = c(1, 1))+
  scale_fill_gradient(low="grey30", high = "grey90")+
  theme_classic(16)+labs(x="observerd", y="fitted")

