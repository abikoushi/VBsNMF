\name{SVBNMF_mtx}
\alias{SVBNMF_mtx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SVBNMF_mtx(file_path, rank,
           b_size = 100, subiter = 1,
           n_epochs = 10, delay = 1,
           forgetting = 0.8, prior_shape = 1, prior_rate = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_path}{
  input data matrix to decompose; path of your \code{.mtx} (matrix market format) file
}
  \item{rank}{
  rank of the output matices
}
  \item{b_size}{
mini-batch size
}
  \item{subiter}{
maximu iteration number in each step
}
  \item{n_epochs}{
number of epochs
}
  \item{delay}{
chooning parameter for learning rate
}
  \item{forgetting}{
chooning parameter for learning rate
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
\seealso{
\code{\link{help}}
}
\examples{
library(Matrix)
library(VBsNMF)

set_data <- function(L, nrow, ncol, center=0, scale=1){
  Z <- matrix(rnorm(L*nrow,0,scale), nrow, L)
  Z <- sweep(Z,1,rowMeans(Z)-center)
  W <- matrix(rnorm(L*ncol,0,scale), L, ncol)
  W <- sweep(W,1,rowMeans(W)-center)
  Y <- matrix(rpois(nrow*ncol, exp(Z) \%*\%exp(W)), nrow, ncol)
  list(Y=Y,Z=Z,W=t(W))
}

set.seed(1234); dat <- set_data(2, 99, 205)

Y <- as(dat$Y, "TsparseMatrix")
writeMM(Y,"test.mtx")

system.time({
out_s <- VBsNMF:::svb_nmf_pois_mtx("test.mtx", rank = 2,
                                   b_size = 10000,
                                   subiter = 5,
                                   forgetting = 0.9, delay = 1,
                                   n_epochs = 100)
})

}