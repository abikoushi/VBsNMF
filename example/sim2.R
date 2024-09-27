library(VBsNMF)
#write file loader

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

