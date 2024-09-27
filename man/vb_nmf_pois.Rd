\name{vb_nmf_pois, svb_nmf_pois}
\alias{vb_nmf_pois}
\alias{svb_nmf_pois}
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
            prior_shape=0.5,
            prior_rate=1e-5)

            
svb_nmf_pois(Y, rank,
            b_size = 100,
            subiter = 10,
            n_epochs = 10,
            delay=1,
            forgetting=0.9,
            prior_shape=0.5,
            prior_rate=0.1)
}
\arguments{
  \item{Y}{
  input matix
}
  \item{rank}{
  rank of the output matices
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Cemgil, A. T. (2009). Bayesian inference for nonnegative matrix factorisation models. Computational intelligence and neuroscience, 2009(1), 785152.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
