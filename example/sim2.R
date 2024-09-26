library(VBsNMF)
set_data <- function(L, nrow, ncol, scale){
  Z <- matrix(rgamma(L*nrow,1,1), nrow, L)
  #Z <- Z/rowSums(Z)
  W <- matrix(rgamma(L*ncol,1,1), L, ncol)
  #W <- W/rowSums(W)
  Y <- matrix(rpois(nrow*ncol, (Z%*%W)*scale), nrow, ncol)
  list(Y=Y,Z=Z,W=t(W))
}

set.seed(1010); dat <- set_data(2, 100, 100, scale = 1)

Y <- as(dat$Y, "TsparseMatrix")
a <- 0.5
b <- 0.1
rank <- 2
alpha_z = matrix(a, Y@Dim[1], rank)
beta_z = matrix(b, 1, rank)
alpha_w = matrix(a, Y@Dim[2], rank)
beta_w = matrix(b, 1, rank)

N1 <- length(Y@x)
rind <-sample.int(N1)

b_size <- 100
sb <- floor(N1/b_size)
m_sb <-b_size*sb
subind <- vector("list",sb+1)
for(i in 1:sb){
  subind[[i]] <- rind[1:b_size+b_size*(i-1)]  
}
subind[[sb+1]] <- rind[(b_size*sb):length(rind)]

learning_rate <- function(t, delay=1, forgetting=0.9){
  (t+delay)^(-forgetting)
}

subiter <- 10
n_epochs <- 10
lp <- numeric(n_epochs)
pb <- txtProgressBar(1, n_epochs, style = 3)
for(ep in 1:n_epochs){
  rho <- learning_rate(ep)
  for(k in 1:length(subind)){
    out_t <- VBsNMF:::doVB_pois_s(Y@x[subind[[k]]],Y@i[subind[[k]]],Y@j[subind[[k]]],
                                  L=2, iter = subiter,
                                  a,b,
                                  N1,
                                  alpha_z, beta_z,
                                  alpha_w, beta_w)
    Ns <- length(subind[[k]])
    alpha_z <- (1-rho)*alpha_z + rho*out_t$shape_row
    beta_z <- (1-rho)*beta_z + rho*out_t$rate_row
    alpha_w <- (1-rho)*alpha_w + rho*out_t$shape_col
    beta_w <- (1-rho)*beta_w + rho*out_t$rate_col
    lp[ep] <- lp[ep] + sum(out_t$logprob[subiter,])
  }  
  setTxtProgressBar(pb,ep)
}

plot(lp, type="l")
Zhat <- sweep(alpha_z,1,c(beta_z),"/")
What <- sweep(alpha_w,1,c(beta_w),"/")
plot(Zhat%*%t(What), dat$Y, col=rgb(0,0,0,0.1))
abline(0,1)


###
a<-0.5
b<-0.1
alpha_z = matrix(a, Y@Dim[1], rank)
beta_z = matrix(b, 1, rank)
alpha_w = matrix(a, Y@Dim[2], rank)
beta_w = matrix(b, 1, rank)
p1 <- length(Y@x)/(100*100)
nrow(alpha_z)
nrow(alpha_w)

out_t <- VBsNMF:::doVB_pois_s(Y@x,Y@i,Y@j,
                              L=2, iter = 10,
                              a,b,
                              p1,
                              alpha_z, beta_z,
                              alpha_w, beta_w)
plot(rowSums(out_t$logprob), type="l")
out$shape_row
out_t$shape_row
out_t$rate_row
Zhat <- with(out_t, sweep(shape_row,1,c(rate_row),"/"))
What <- with(out_t, sweep(shape_col,1,c(rate_col),"/"))
plot(Zhat%*%t(What), Y, col=rgb(0,0,0,0.1))
abline(0,1)

system.time({
  out <- vb_nmf_pois(Y, rank=2, iter=10, prior_rate=0.1)
})

Zhat <- with(out, sweep(shape_row,1,c(rate_row),"/"))
What <- with(out, sweep(shape_col,1,c(rate_col),"/"))
plot(Zhat%*%t(What), Y, col=rgb(0,0,0,0.1))
abline(0,1)




out_s <- svb_nmf_pois(Y, rank=2)
plot(out_s$elbo)
class(out_t)
