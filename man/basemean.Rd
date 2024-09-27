\name{basemean, coefmean}?bas
\alias{basemean}
\alias{coefmean}
\title{
The posterior mean of base or coefficient matrix
}
\description{
\code{basemean} returns the posterior mean of base matrix.
\code{coefmean} returns the posterior mean of coefficient matrix.
}
\usage{
basemean(obj)

coefmean(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  output object of \code{vb_nmf_pois} or \code{svb_nmf_pois}
}
}
\details{

}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (obj) 
{
    stopifnot(any(class(obj) == "nmf_pois_posterior"))
    with(obj, sweep(shape_row, 1, c(rate_row), "/"))
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
