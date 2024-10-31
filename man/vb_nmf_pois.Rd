\name{vb_nmf_pois}
\alias{vb_nmf_pois}
\alias{em_nmf_pois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variational Bayes method for sparase nonegative matrix facttorization
}
\description{
Efficient implementation enjoying the sparsity of the Variational Bayes method for nonnegative matrix factorization.
}
\usage{
vb_nmf_pois(Y, rank,
            iter=100,
            prior_shape=1,
            prior_rate=1)
            
em_nmf_pois(Y, rank,
            iter=100,
            prior_shape=1,
            prior_rate=1)
}
\arguments{
  \item{Y}{
  input data matrix to decompose
}
  \item{rank}{
  rank of the output matices
}
  \item{iter}{
  the number of the iterations
}
  \item{prior_shape}{
  shape paremeter of the gamma prior distribution
}
  \item{prior_rate}{
  rate paremeter of the gamma prior distribution
}
}
\value{
  \item{shape_row }{shape parameters of the posterir gamma distribution for row-axis}
  \item{rate_row }{rate parameters of the posterir gamma distribution for row-axis}
  \item{shape_col }{shape parameters of the posterir gamma distribution for col-axis}
  \item{rate_col }{rate parameters of the posterir gamma distribution for col-axis}
  \item{logprob }{lower bound of the log-joint probability}
}
\references{
Cemgil, A. T. (2009). Bayesian inference for nonnegative matrix factorisation models. Computational intelligence and neuroscience, 2009(1), 785152.
}
\examples{
#draw random matrix
set_data <- function(L, nrow, ncol, center=0, scale=1){
  Z <- matrix(rnorm(L*nrow,0,scale), nrow, L)
  Z <- sweep(Z,1,rowMeans(Z)-center)
  W <- matrix(rnorm(L*ncol,0,scale), L, ncol)
  W <- sweep(W,1,rowMeans(W)-center)
  Y <- matrix(rpois(nrow*ncol, exp(Z)%*%exp(W)), nrow, ncol)
  list(Y=Y,Z=Z,W=t(W))
}

set.seed(1234); dat <- set_data(2, 99, 205)
indna <- sample.int(length(dat$Y), size = 100)
Y <- dat$Y
Y[indna] <- NA

Y <- as(Y, "TsparseMatrix")
hist(dat$Y, breaks = "FD")

system.time({
  out <- vb_nmf_pois(dat$Y, rank=2, iter=100, prior_rate=1)
})

plot(out$logprob)

Zhat <- basemean(out)
What <- coefmean(out)

fit <- Zhat\%*\%What

basecol <-rgb(0,0,0,0.1)
plot(dat$Y[-indna], fit[-indna], col=basecol, cex=0.5)
points(dat$Y[indna], fit[indna], col=rgb(1,0.5,0,0.3), pch=16)
abline(0,1, col="lightgrey", lty=2)
}
