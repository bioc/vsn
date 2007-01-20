\name{logLik-methods}
\alias{logLik-methods}
\alias{logLik,vsnInput-method}

\title{Calculate the log likelihood and its gradient for
  the vsn model}
\description{Calculate the log likelihood and its gradient
  for the vsn model}
\usage{\S4method{logLik}{vsnInput}(object, p, refh, refsigma)}
\arguments{
  \item{object}{A \code{\link{vsnInput}} object.}
  \item{p}{A matrix whose columns are the set of parameters at which
    the likelihoods are to be evaluated.}
  \item{refh}{Numeric vector of the same length as \code{nrow(object)}.
    The arguments \code{refh} and \code{refsigma} must be either both present or
    both missing.}
  \item{refsigma}{Numeric scalar.}
}
\value{
  A numeric matrix of size \code{nrow(p)+1} by \code{ncol(p)}.
  Its columns correspond to the columns of \code{p}.
  Its first row are the likelihood values, its rows \code{2...nrow(p)+1}
  contain the gradients.

  If \code{refh} and \code{refsigma} are
  specified, the ordinary negative log likelihood is calculated using these
  parameters as given. If they are not specified, the profile negative log likelihood
  is calculated.
}

\details{
  This is an R interface to the likelihood computations in vsn (which are done in C).
  It is mostly useful for testing purposes.
}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}
\seealso{\code{\link{vsn2}}}

\examples{

data("kidney")

v = new("vsnInput", x=exprs(kidney),
  pstart=array(as.numeric(NA), dim=c(1, ncol(kidney), 2)))
 
fit = vsn2(kidney)

p = sapply(exp(seq(-1, 1, length=31)), function(f) fit@par*c(1,1,f,1))

ll = logLik(v, p)

plot(p[3, ], ll[1, ], type="l", xlab=expression(b[1]), ylab=expression(-log(L)), log="x")
abline(v=fit@par[3], col="red")
}

\keyword{}
