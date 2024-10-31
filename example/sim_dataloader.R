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

set.seed(1234); dat <- set_data(2, 99, 205)
mean(dat$Y==0)
Y <- as(dat$Y, "TsparseMatrix")
writeMM(Y,"test.mtx")

system.time({
out_s <- VBsNMF:::svb_nmf_pois_mtx("test.mtx", rank = 2,
                                   b_size = 10000,
                                   subiter = 5,
                                   forgetting = 0.9, delay = 1,
                                   n_epochs = 100)
})
hist(dat$Y, breaks = "FD")
plot(out_s$logprob, type="l")
Zhat <- basemean(out_s)
What <- coefmean(out_s)

fit <- Zhat%*%What
ggplot(data=NULL, aes(x=c(dat$Y), y=c(fit)))+
  geom_abline(intercept = 0, slope = 1, col="lightgrey", linetype=2) +
  geom_bin2d(aes(fill=after_stat(log10(count))))+
  theme_classic(16)+labs(x="observerd", y="fitted")

