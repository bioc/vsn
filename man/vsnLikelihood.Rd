\name{logLik-methods}
\alias{logLik-methods}
\alias{logLik,vsnInput-method}
\alias{plotVsnLogLik}

\title{Calculate the log likelihood and its gradient for
  the vsn model}
\description{\code{logLik} calculates the log likelihood and its gradient
  for the vsn model. \code{plotVsnLogLik} makes a false color plot for
  a 2D section of the likelihood landscape.}
\usage{
\S4method{logLik}{vsnInput}(object, p, refh=numeric(0), refsigma=as.numeric(NA))

plotVsnLogLik(object, p, whichp=1:2, expand=1, ngrid=31, fun=logLik, ...)
}
\arguments{
  \item{object}{A \code{\link{vsnInput}} object.}
  \item{p}{For \code{plotVsnLogLik}, a vector or a 3D array with the
    point in parameter space around which to plot the likelihood.
    For \code{logLik}, a matrix whose columns are the set of parameters at which
    the likelihoods are to be evaluated.}
  \item{refh}{Numeric vector of length 0 or \code{nrow(object)}.
    If the length is 0, there is no reference and and \code{refsigma} must be
    missing. See \code{\link{vsn2}}.}
  \item{refsigma}{Numeric scalar.}
  \item{whichp}{Numeric vector of length 2, with the indices of those
    two parameters in \code{p} along which the section is to be taken.}
  \item{expand}{Numeric vector of length 1 or 2 with expansion factors
    for the plot range. The range is auto-calculated using a heuristic,
    but that heuristic is not very good in particular for the offset (a)
    parameters; see example.}
  \item{ngrid}{Grid size.}
  \item{fun}{Function to use for log-likelihood calculation. This
    parameter is exposed only for testing purposes.}
  \item{...}{Arguments that get passed on to \code{fun}, use this for
    \code{refh} and \code{refsigma}.}
}
\value{
For \code{logLik}, a numeric matrix of size \code{nrow(p)+1} by \code{ncol(p)}.
Its columns correspond to the columns of \code{p}.
Its first row are the likelihood values, its rows \code{2...nrow(p)+1}
contain the gradients.
If \code{refh} and \code{refsigma} are
specified, the ordinary negative log likelihood is calculated using these
parameters as given. If they are not specified, the profile negative log likelihood
is calculated.

For \code{plotVsnLogLik}, a dataframe with the 2D grid coordinates and
log likelihood values.
}

\details{
  \code{logLik} is an R interface to the likelihood computations in vsn (which are done in C).
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
abline(v=coef(fit)[3], col="red")

plotVsnLogLik(v, coef(fit), whichp=c(1,3), expand=c(0.01, 1))
}

\keyword{}
