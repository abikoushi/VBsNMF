library(Matrix)
library(ggplot2)
library(VBsNMF)

#par <- "VBparams/shape_col.csv"
par <- "test.mtx"
size <- VBsNMF:::size_mtx(par)
out_t <- VBsNMF:::doVB_pois_s_mtx(par, 2, iter = 100, subiter = 5,
                         a = 1, b = 1, N1 = size[3],
                         Nr = size[1], Nc = size[2],
                         ns = 1000)

plot(out_t$logprob, type = "l")
fit <-  basemean(out_t)%*%coefmean(out_t)

Y <- readMM(par)
plot(as.matrix(Y), fit, pch=".")
abline(0,1,col="cornflowerblue")

set_data <- function(L, nrow, ncol, center=0, scale=1){
  Z <- matrix(rnorm(L*nrow,0,scale), nrow, L)
  Z <- sweep(Z,1,rowMeans(Z)-center)
  W <- matrix(rnorm(L*ncol,0,scale), L, ncol)
  W <- sweep(W,1,rowMeans(W)-center)
  Y <- matrix(rpois(nrow*ncol, exp(Z)%*%exp(W)), nrow, ncol)
  list(Y=Y,Z=Z,W=t(W))
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

system.time({
out_s <- VBsNMF:::SVBNMF("test.mtx",
                         rank = 2,
                         b_size = 100,
                         subiter = 5,
                         forgetting = 0.9, delay = 1,
                         n_epochs = 200,
                         setinit = VBsNMF:::setinit_csv,
                         updatefun = VBsNMF:::updatefun_csv, 
                         draw_param = VBsNMF:::draw_param_csv)
})

plot(out_s$logprob, type="l")

out_s2 <- list(
  shape_col = as.matrix(read.csv("VBparams/shape_col.csv", header = FALSE)),
  shape_row = as.matrix(read.csv("VBparams/shape_row.csv", header = FALSE)),
  rate_col = as.matrix(read.csv("VBparams/rate_col.csv", header = FALSE)),
  rate_row = as.matrix(read.csv("VBparams/rate_row.csv", header = FALSE))
  )

Zhat <- basemean(out_s2)
What <- coefmean(out_s2)
fit_s <- Zhat%*%What

colSums(Zhat)
rowSums(What)

# Y <- readMM("test.mtx")
# plot(fit_s, as.matrix(Y), pch=".")
ggplot(data=NULL, aes(x=c(dat$Y)))+
  geom_abline(intercept = 0, slope = 1, col="grey", linetype=1) +
  geom_jitter(aes(y=c(fit), colour = "VB"), size=0.1, alpha=0.1, width = 0.5, height = 0.5)+
  geom_jitter(aes(y=c(fit_s), colour = "SVB"), size=0.1, alpha=0.1, width = 0.5, height = 0.5)+
  guides(colour=guide_legend(override.aes = list(alpha=1, size=1)))+
  theme_classic(16)+labs(x="observerd", y="fitted", colour="method")

ggplot(data=NULL, aes(x=c(dat$Y), y=c(fit_s)))+
  geom_abline(intercept = 0, slope = 1, col="grey", linetype=1) +
  geom_bin2d(aes(fill = after_stat(log2(count))))+
  theme_classic(16)+labs(x="observerd", y="fitted")

