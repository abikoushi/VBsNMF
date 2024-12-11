library(Matrix)
library(ggplot2)
library(VBsNMF)
library(rbenchmark)
library(pryr)
library(bench)

size <- VBsNMF:::size_mtx("test.mtx")
# bag <- sort(sample.int(size[3], 100))
# 
# #VBsNMF:::readmtx("test.mtx",bag)
# benchmark(res1 = VBsNMF:::dataloader_mtx(file_path = "test.mtx", bag = bag),
#           res2 = VBsNMF:::read_mtx(readtxt = "test.mtx", bag = bag))
# res1 = VBsNMF:::dataloader_mtx(file_path = "test.mtx", bag = bag)
# res2 = VBsNMF:::read_mtx(readtxt = "test.mtx", bag = bag)
# res1[[1]]
# c(res2[[1]])

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
bench::mark({
out_s <- VBsNMF:::SVBNMF("test.mtx",
                         rank = 2,
                         b_size = 10000,
                         subiter = 5,
                         forgetting = 0.8, delay = 1,
                         n_epochs = 200)
})

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

# ca005 <- rgb(0,0,0,0.01)
# Y <- readMM("test.mtx")
# plot(fit_s, as.matrix(Y), col=ca005, cex=0.5)
# abline(0,1, col="grey")
